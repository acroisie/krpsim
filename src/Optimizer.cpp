#include "Optimizer.hpp"
#include <algorithm>
#include <random>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_set>

Optimizer::Optimizer(const Config& config, int timeLimit)
    : config_(config), timeLimit_(timeLimit), simulator_(config) {
}

Optimizer::Solution Optimizer::optimize() {
    std::vector<Solution> candidates;
    
    // 1. Essai de diverses stratégies
    candidates.push_back(greedyOptimize());
    candidates.push_back(randomSearch(100));
    
    // Choix d'une stratégie initiale pour démarrer la recherche locale
    std::sort(candidates.begin(), candidates.end(), compareSolutions);
    
    Solution bestSolution = candidates.front();
    std::cout << "Initial best score: " << bestSolution.score << std::endl;
    
    // 2. Améliorations par recherche locale
    bestSolution = simulatedAnnealing(bestSolution, 1000);
    std::cout << "After simulated annealing: " << bestSolution.score << std::endl;
    
    bestSolution = hillClimbing(bestSolution, 500);
    std::cout << "After hill climbing: " << bestSolution.score << std::endl;
    
    return bestSolution;
}

void Optimizer::printResults(const Solution& solution) {
    Simulator::Result result = simulator_.simulate(solution.processSequence, timeLimit_);
    
    std::cout << "Nice file! " << config_.getProcesses().size() << " processes, "
              << config_.getStocks().size() << " stocks, "
              << config_.getOptimizeGoal().size() << " to optimize" << std::endl;
    
    std::cout << "Main walk:" << std::endl;
    for (const auto& log : result.executionLog) {
        std::cout << log.first << ":" << log.second << std::endl;
    }
    
    std::cout << "no more process doable at time " << result.finalTime << std::endl;
    std::cout << "Stock :" << std::endl;
    for (const auto& stock : result.finalStocks) {
        std::cout << stock.first << " => " << stock.second << std::endl;
    }
}

Optimizer::Solution Optimizer::greedyOptimize() {
    // Stratégie 1: Priorité basée sur l'efficacité des ressources
    Solution solution;
    solution.processSequence = prioritizeByResourceEfficiency();
    solution.score = evaluate(solution);
    
    // Stratégie 2: Analyse des dépendances
    Solution depSolution;
    depSolution.processSequence = analyzeProcessDependencies();
    depSolution.score = evaluate(depSolution);
    
    return (solution.score > depSolution.score) ? solution : depSolution;
}

std::vector<std::string> Optimizer::analyzeProcessDependencies() {
    // Graphe de dépendances entre processus
    std::map<std::string, std::vector<std::string>> dependsOn;
    std::map<std::string, std::vector<std::string>> produces;
    std::map<std::string, int> producedBy;
    
    // Identifier les ressources produites par chaque processus
    for (const auto& process : config_.getProcesses()) {
        for (const auto& output : process.outputs) {
            produces[process.name].push_back(output.first);
            producedBy[output.first]++;
        }
    }
    
    // Identifier les dépendances entre processus
    for (const auto& process : config_.getProcesses()) {
        for (const auto& input : process.inputs) {
            // Si l'input est produit par un autre processus
            if (producedBy.find(input.first) != producedBy.end()) {
                for (const auto& otherProcess : config_.getProcesses()) {
                    if (otherProcess.name == process.name) continue;
                    
                    for (const auto& output : otherProcess.outputs) {
                        if (output.first == input.first) {
                            dependsOn[process.name].push_back(otherProcess.name);
                        }
                    }
                }
            }
        }
    }
    
    // Classer les processus en fonction des dépendances
    std::vector<std::string> result;
    std::unordered_set<std::string> visited;
    
    // Ajouter d'abord les processus sans dépendances
    for (const auto& process : config_.getProcesses()) {
        if (dependsOn.find(process.name) == dependsOn.end() || dependsOn[process.name].empty()) {
            result.push_back(process.name);
            visited.insert(process.name);
        }
    }
    
    // Puis les autres en respectant les dépendances
    while (visited.size() < config_.getProcesses().size()) {
        for (const auto& process : config_.getProcesses()) {
            if (visited.find(process.name) != visited.end()) continue;
            
            bool allDependenciesSatisfied = true;
            if (dependsOn.find(process.name) != dependsOn.end()) {
                for (const auto& dep : dependsOn[process.name]) {
                    if (visited.find(dep) == visited.end()) {
                        allDependenciesSatisfied = false;
                        break;
                    }
                }
            }
            
            if (allDependenciesSatisfied) {
                result.push_back(process.name);
                visited.insert(process.name);
            }
        }
        
        // Si aucun processus n'a été ajouté, briser un cycle
        if (visited.size() < config_.getProcesses().size() && 
            visited.size() == result.size()) {
            for (const auto& process : config_.getProcesses()) {
                if (visited.find(process.name) == visited.end()) {
                    result.push_back(process.name);
                    visited.insert(process.name);
                    break;
                }
            }
        }
    }
    
    return result;
}

std::vector<std::string> Optimizer::prioritizeByResourceEfficiency() {
    std::vector<std::pair<std::string, double>> processEfficiency;
    
    // Déterminer les objectifs d'optimisation
    const auto& goals = config_.getOptimizeGoal();
    bool optimizeTime = goals.empty() || (goals.size() == 1 && goals[0] == "time");
    
    for (const auto& process : config_.getProcesses()) {
        double efficiency = 0.0;
        
        // Valeur des outputs
        double outputValue = 0.0;
        for (const auto& output : process.outputs) {
            if (!optimizeTime && std::find(goals.begin(), goals.end(), output.first) != goals.end()) {
                // Si c'est une ressource à optimiser, elle a une valeur plus élevée
                outputValue += output.second * 10.0;
            } else {
                outputValue += output.second;
            }
        }
        
        // Coût des inputs
        double inputCost = 0.0;
        for (const auto& input : process.inputs) {
            inputCost += input.second;
        }
        
        // Temps d'exécution
        double timeWeight = optimizeTime ? 5.0 : 1.0;
        
        // Efficacité = (valeur des outputs) / (coût des inputs * temps)
        if (inputCost > 0 && process.nbCycle > 0) {
            efficiency = outputValue / (inputCost * process.nbCycle / timeWeight);
        } else if (process.nbCycle > 0) {
            efficiency = outputValue / (process.nbCycle / timeWeight);
        } else {
            efficiency = outputValue * 1000; // Très efficace si temps = 0
        }
        
        processEfficiency.push_back({process.name, efficiency});
    }
    
    // Trier par efficacité décroissante
    std::sort(processEfficiency.begin(), processEfficiency.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Extraire uniquement les noms de processus
    std::vector<std::string> result;
    for (const auto& p : processEfficiency) {
        result.push_back(p.first);
    }
    
    return result;
}

double Optimizer::evaluate(const Solution& solution) {
    Simulator::Result result = simulator_.simulate(solution.processSequence, timeLimit_, false);
    return result.score;
}

Optimizer::Solution Optimizer::randomSolution() {
    Solution solution;
    
    // Créer une permutation aléatoire des processus
    for (const auto& process : config_.getProcesses()) {
        solution.processSequence.push_back(process.name);
    }
    
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(solution.processSequence.begin(), solution.processSequence.end(), g);
    
    solution.score = evaluate(solution);
    return solution;
}

Optimizer::Solution Optimizer::randomSearch(int numTrials) {
    Solution bestSolution;
    bestSolution.score = std::numeric_limits<double>::lowest();
    
    for (int i = 0; i < numTrials; ++i) {
        Solution candidate = randomSolution();
        if (candidate.score > bestSolution.score) {
            bestSolution = candidate;
        }
    }
    
    return bestSolution;
}

Optimizer::Solution Optimizer::neighborSolution(const Solution& solution) {
    Solution neighbor = solution;
    
    // Générer une solution voisine en permutant deux processus aléatoires
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, solution.processSequence.size() - 1);
    
    int i = dist(gen);
    int j = dist(gen);
    while (i == j && solution.processSequence.size() > 1) {
        j = dist(gen);
    }
    
    std::swap(neighbor.processSequence[i], neighbor.processSequence[j]);
    
    neighbor.score = evaluate(neighbor);
    return neighbor;
}

bool Optimizer::compareSolutions(const Solution& a, const Solution& b) {
    return a.score > b.score;
}

Optimizer::Solution Optimizer::hillClimbing(const Solution& initialSolution, int iterations) {
    Solution current = initialSolution;
    
    for (int i = 0; i < iterations; ++i) {
        Solution neighbor = neighborSolution(current);
        
        if (neighbor.score > current.score) {
            current = neighbor;
        }
    }
    
    return current;
}

Optimizer::Solution Optimizer::simulatedAnnealing(const Solution& initialSolution, int iterations) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);
    
    Solution current = initialSolution;
    Solution best = current;
    
    // Paramètres de température
    double temperature = 1.0;
    double coolingRate = 0.99;
    
    for (int i = 0; i < iterations; ++i) {
        Solution neighbor = neighborSolution(current);
        
        // Décider si on accepte la nouvelle solution
        if (neighbor.score > current.score) {
            current = neighbor;
            
            if (current.score > best.score) {
                best = current;
            }
        } else {
            // Accepter des solutions moins bonnes avec une probabilité qui diminue avec la température
            double acceptanceProbability = exp((neighbor.score - current.score) / temperature);
            if (dist(gen) < acceptanceProbability) {
                current = neighbor;
            }
        }
        
        // Refroidissement
        temperature *= coolingRate;
    }
    
    return best;
}