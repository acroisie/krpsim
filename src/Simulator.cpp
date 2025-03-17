#include "Simulator.hpp"
#include <algorithm>
#include <limits>
#include <iostream>

Simulator::Simulator(const Config& config) 
    : config_(config), currentTime_(0) {
    // Initialisation des stocks
    for (const auto& stock : config.getStocks()) {
        stocks_[stock.name] = stock.quantity;
    }
}

Simulator::Result Simulator::simulate(
    const std::vector<std::string>& processSequence, 
    int timeLimit,
    bool generateLog) {
    
    // Réinitialisation de l'état
    stocks_.clear();
    for (const auto& stock : config_.getStocks()) {
        stocks_[stock.name] = stock.quantity;
    }
    eventQueue_ = std::priority_queue<Event, std::vector<Event>, std::greater<Event>>();
    executionLog_.clear();
    currentTime_ = 0;
    
    // Mapping des priorités des processus
    std::map<std::string, int> processPriority;
    for (size_t i = 0; i < processSequence.size(); ++i) {
        processPriority[processSequence[i]] = i;
    }
    
    bool processStarted = true;
    while (currentTime_ < timeLimit && processStarted) {
        processStarted = false;
        
        // Traitement des événements de fin de processus
        while (!eventQueue_.empty() && eventQueue_.top().time <= currentTime_) {
            Event event = eventQueue_.top();
            eventQueue_.pop();
            
            if (!event.isStart) { // Événement de fin de processus
                const Process* process = getProcessByName(event.processName);
                if (process) {
                    produceResources(*process);
                }
            }
        }
        
        // Recherche du processus exécutable avec la priorité la plus élevée
        const Process* bestProcess = nullptr;
        int bestPriority = std::numeric_limits<int>::max();
        
        for (const auto& process : config_.getProcesses()) {
            if (canRunProcess(process)) {
                int priority = processPriority[process.name];
                if (priority < bestPriority) {
                    bestPriority = priority;
                    bestProcess = &process;
                }
            }
        }
        
        // Lancement du processus sélectionné
        if (bestProcess) {
            if (generateLog) {
                executionLog_.push_back({currentTime_, bestProcess->name});
            }
            startProcess(*bestProcess, currentTime_);
            processStarted = true;
        }
        
        // Avancement du temps
        if (!processStarted && !eventQueue_.empty()) {
            // Si aucun processus n'a démarré mais qu'il en reste en cours d'exécution,
            // on avance jusqu'au prochain événement
            currentTime_ = eventQueue_.top().time;
        } else {
            // Sinon on avance d'une unité de temps
            currentTime_++;
        }
    }
    
    // Attendre la fin de tous les processus en cours
    while (!eventQueue_.empty()) {
        Event event = eventQueue_.top();
        eventQueue_.pop();
        
        if (event.time > timeLimit) {
            break;  // On s'arrête si on dépasse la limite de temps
        }
        
        currentTime_ = event.time;
        
        if (!event.isStart) { // Événement de fin de processus
            const Process* process = getProcessByName(event.processName);
            if (process) {
                produceResources(*process);
            }
        }
    }
    
    // Construction du résultat
    Result result;
    result.executionLog = executionLog_;
    result.finalStocks = stocks_;
    result.finalTime = currentTime_;
    result.score = calculateScore(stocks_, currentTime_);
    
    return result;
}

double Simulator::calculateScore(const std::map<std::string, int>& stocks, int time) const {
    const auto& goals = config_.getOptimizeGoal();
    
    // Par défaut, on minimise le temps
    if (goals.empty() || (goals.size() == 1 && goals[0] == "time")) {
        return -time;
    }
    
    // Sinon, on maximise la ressource indiquée
    double score = 0.0;
    for (const auto& goal : goals) {
        if (goal != "time") {
            auto it = stocks.find(goal);
            if (it != stocks.end()) {
                score += it->second;
            }
        }
    }
    
    return score;
}

bool Simulator::canRunProcess(const Process& process) const {
    for (const auto& input : process.inputs) {
        auto it = stocks_.find(input.first);
        if (it == stocks_.end() || it->second < input.second) {
            return false;
        }
    }
    return true;
}

void Simulator::startProcess(const Process& process, int startTime) {
    // Consommer les ressources nécessaires
    consumeResources(process);
    
    // Ajouter l'événement de fin de processus à la file d'attente
    Event endEvent;
    endEvent.time = startTime + process.nbCycle;
    endEvent.processName = process.name;
    endEvent.isStart = false;
    eventQueue_.push(endEvent);
}

void Simulator::consumeResources(const Process& process) {
    for (const auto& input : process.inputs) {
        stocks_[input.first] -= input.second;
    }
}

void Simulator::produceResources(const Process& process) {
    for (const auto& output : process.outputs) {
        stocks_[output.first] += output.second;
    }
}

const Process* Simulator::getProcessByName(const std::string& name) const {
    for (const auto& process : config_.getProcesses()) {
        if (process.name == name) {
            return &process;
        }
    }
    return nullptr;
}