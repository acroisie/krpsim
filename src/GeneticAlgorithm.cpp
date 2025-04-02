#include "GeneticAlgorithm.hpp"
#include <vector>
#include <string>
#include <random>
#include <algorithm> // Pour std::sort, std::min, std::max, std::max_element, std::all_of
#include <iostream>  // Pour l'affichage (debug/info)
#include <limits>    // Pour std::numeric_limits
#include <cmath>     // Pour std::abs, std::isfinite

// Constructeur (inchangé)
GeneticAlgorithm::GeneticAlgorithm(const Config& config,
                                   Simulator& simulator,
                                   int populationSize,
                                   double mutationRate,
                                   double crossoverRate,
                                   int eliteCount,
                                   int minSequenceLength,
                                   int maxSequenceLength)
    : systemConfig_(config),
      planningSimulator_(simulator),
      populationSize_(populationSize),
      mutationRate_(mutationRate),
      crossoverRate_(crossoverRate),
      eliteCount_(std::max(0, std::min(eliteCount, populationSize))),
      minSequenceLength_(minSequenceLength),
      maxSequenceLength_(maxSequenceLength),
      randomGenerator_(std::random_device{}())
{
    availableProcessNames_.clear();
    for (const auto& proc : systemConfig_.getProcesses()) {
        availableProcessNames_.push_back(proc.name);
    }
    if (availableProcessNames_.empty()) {
        throw std::runtime_error("GeneticAlgorithm Error: No processes defined in the configuration.");
    }
     if (minSequenceLength_ > maxSequenceLength_ || minSequenceLength_ <= 0) {
         throw std::runtime_error("GeneticAlgorithm Error: Invalid sequence length parameters.");
     }
}

// createRandomIndividual (inchangé)
Individual GeneticAlgorithm::createRandomIndividual() {
    std::uniform_int_distribution<> lenDist(minSequenceLength_, maxSequenceLength_);
    std::uniform_int_distribution<> procIdxDist(0, static_cast<int>(availableProcessNames_.size()) - 1);
    int sequenceLength = lenDist(randomGenerator_);
    std::vector<std::string> sequence;
    sequence.reserve(sequenceLength);
    for (int j = 0; j < sequenceLength; ++j) {
        sequence.push_back(availableProcessNames_[procIdxDist(randomGenerator_)]);
    }
    return Individual(sequence);
}

// initializePopulation (inchangé)
void GeneticAlgorithm::initializePopulation() {
    currentPopulation_.clear();
    currentPopulation_.reserve(populationSize_);
    for (int i = 0; i < populationSize_; ++i) {
        currentPopulation_.push_back(createRandomIndividual());
    }
    std::cout << "Population initialized with " << populationSize_ << " individuals." << std::endl;
}

// Évalue le fitness de tous les individus de la population actuelle (Non-const)
void GeneticAlgorithm::evaluatePopulationFitness() {
    for (auto& individual : currentPopulation_) {
        // *** APPEL CORRIGÉ ICI ***
        // Appeler la nouvelle fonction qui retourne SimulationResult
        Simulator::SimulationResult result = planningSimulator_.calculateFitnessAndLogs(
            individual.processExecutionAttemptSequence);
        // Extraire le score de fitness du résultat
        individual.fitnessScore = result.fitness;

        // Debug: Afficher les scores calculés
        // std::cout << "  Evaluated fitness: " << individual.fitnessScore << std::endl;
    }
     // Debug: Vérifier si tous les scores sont identiques et très bas
     /*
     if (!currentPopulation_.empty()) {
         double firstScore = currentPopulation_[0].fitnessScore;
         bool allSame = std::all_of(currentPopulation_.begin(), currentPopulation_.end(),
                                    [firstScore](const Individual& ind){ return ind.fitnessScore == firstScore; });
         if (allSame && firstScore < -1e17) {
             std::cout << "Debug: All individuals have the same minimum fitness score." << std::endl;
         }
     }
     */
}

// Sélectionne les parents par la méthode de la roue de la fortune (Non-const)
std::vector<size_t> GeneticAlgorithm::selectParentsViaRouletteWheel() {
     std::vector<size_t> selectedIndices;
     if (currentPopulation_.empty()) {
         return selectedIndices;
     }

     double minFitness = std::numeric_limits<double>::max();
     double maxFitness = std::numeric_limits<double>::lowest();

     for (const auto& indiv : currentPopulation_) {
         // Ignorer les fitness invalides pour min/max
         if (std::isfinite(indiv.fitnessScore)) {
             if (indiv.fitnessScore < minFitness) minFitness = indiv.fitnessScore;
             if (indiv.fitnessScore > maxFitness) maxFitness = indiv.fitnessScore;
         }
     }

     // Gérer le cas où aucun fitness valide n'a été trouvé
     if (minFitness == std::numeric_limits<double>::max() || maxFitness == std::numeric_limits<double>::lowest()) {
         std::cerr << "Warning: No valid fitness scores found in population. Selecting parents randomly." << std::endl;
         selectedIndices.reserve(populationSize_);
         std::uniform_int_distribution<size_t> randDist(0, currentPopulation_.empty() ? 0 : currentPopulation_.size() - 1);
         for(int i = 0; i < populationSize_; ++i) {
             if (currentPopulation_.empty()) break;
             selectedIndices.push_back(randDist(randomGenerator_));
         }
         return selectedIndices;
     }


     // Si tous les fitness valides sont identiques
     if (minFitness == maxFitness) {
          std::cout << "Warning: All individuals have the same valid fitness (" << minFitness << "). Selecting parents randomly." << std::endl;
          selectedIndices.reserve(populationSize_);
          std::uniform_int_distribution<size_t> randDist(0, currentPopulation_.size() - 1);
          for(int i = 0; i < populationSize_; ++i) {
              selectedIndices.push_back(randDist(randomGenerator_));
          }
          return selectedIndices;
     }

     // Normalisation des scores entre 0.001 et ~1.001
     double range = maxFitness - minFitness;
     if (range <= 0 || !std::isfinite(range)) { range = 1.0; } // Sécurité

     std::vector<double> normalizedFitnessValues;
     normalizedFitnessValues.reserve(currentPopulation_.size());
     double totalNormalizedFitness = 0.0;

     for (const auto& indiv : currentPopulation_) {
         double normalized = 0.001; // Score minimum si fitness invalide ou égal au min
         if (std::isfinite(indiv.fitnessScore)) {
              normalized = (indiv.fitnessScore - minFitness) / range + 0.001;
         }
         normalizedFitnessValues.push_back(normalized);
         totalNormalizedFitness += normalized;
     }


     if (totalNormalizedFitness <= 0.0 || !std::isfinite(totalNormalizedFitness)) {
         std::cerr << "Warning: Total normalized fitness is non-positive or invalid (" << totalNormalizedFitness << "). Selecting parents randomly." << std::endl;
          selectedIndices.reserve(populationSize_);
          std::uniform_int_distribution<size_t> randDist(0, currentPopulation_.empty() ? 0 : currentPopulation_.size() - 1);
          for(int i = 0; i < populationSize_; ++i) {
               if (currentPopulation_.empty()) break;
              selectedIndices.push_back(randDist(randomGenerator_));
          }
          return selectedIndices;
     }

     // Distribution et sélection par roue
     std::uniform_real_distribution<> dist(0.0, totalNormalizedFitness);
     selectedIndices.reserve(populationSize_);
     for (int i = 0; i < populationSize_; ++i) {
         double randomPoint = dist(randomGenerator_);
         double currentSum = 0.0;
         bool selected = false;
         for (size_t j = 0; j < currentPopulation_.size(); ++j) {
             currentSum += normalizedFitnessValues[j];
             if (randomPoint <= currentSum) {
                 selectedIndices.push_back(j);
                 selected = true;
                 break;
             }
         }
         if (!selected && !currentPopulation_.empty()) {
             selectedIndices.push_back(currentPopulation_.size() - 1);
         } else if (!selected && currentPopulation_.empty()){
              break;
         }
     }
     return selectedIndices;
}


// performCrossover (inchangé)
std::pair<Individual, Individual> GeneticAlgorithm::performCrossover(const Individual& parent1, const Individual& parent2) {
    std::uniform_real_distribution<> crossDist(0.0, 1.0);
    if (crossDist(randomGenerator_) > crossoverRate_) {
        return {parent1, parent2};
    }
    const auto& seq1 = parent1.processExecutionAttemptSequence;
    const auto& seq2 = parent2.processExecutionAttemptSequence;
    size_t len1 = seq1.size();
    size_t len2 = seq2.size();
    if (len1 < 2 || len2 < 2) { return {parent1, parent2}; }
    size_t shorterLength = std::min(len1, len2);
    std::uniform_int_distribution<size_t> pointDist(0, shorterLength - 1);
    size_t point1 = pointDist(randomGenerator_);
    size_t point2 = pointDist(randomGenerator_);
    if (point1 == point2) { point2 = (point1 + 1) % shorterLength; }
    if (point1 > point2) { std::swap(point1, point2); }
    std::vector<std::string> childSeq1, childSeq2;
    childSeq1.reserve(len1); childSeq2.reserve(len2);
    childSeq1.insert(childSeq1.end(), seq1.begin(), seq1.begin() + point1);
    if (point2 < len2) { childSeq1.insert(childSeq1.end(), seq2.begin() + point1, seq2.begin() + point2); }
    else if (point1 < len2) { childSeq1.insert(childSeq1.end(), seq2.begin() + point1, seq2.end()); }
    if (point2 < len1) { childSeq1.insert(childSeq1.end(), seq1.begin() + point2, seq1.end()); }
    childSeq2.insert(childSeq2.end(), seq2.begin(), seq2.begin() + point1);
    if (point2 < len1) { childSeq2.insert(childSeq2.end(), seq1.begin() + point1, seq1.begin() + point2); }
    else if (point1 < len1) { childSeq2.insert(childSeq2.end(), seq1.begin() + point1, seq1.end()); }
    if (point2 < len2) { childSeq2.insert(childSeq2.end(), seq2.begin() + point2, seq2.end()); }
    return {Individual(childSeq1), Individual(childSeq2)};
}

// applyMutation (inchangé)
Individual GeneticAlgorithm::applyMutation(const Individual& individual) {
    Individual mutatedIndividual = individual;
    std::uniform_real_distribution<> probDist(0.0, 1.0);
    std::uniform_int_distribution<> procIdxDist(0, static_cast<int>(availableProcessNames_.size()) - 1);
    for (size_t i = 0; i < mutatedIndividual.processExecutionAttemptSequence.size(); ++i) {
        if (probDist(randomGenerator_) < mutationRate_) {
            mutatedIndividual.processExecutionAttemptSequence[i] = availableProcessNames_[procIdxDist(randomGenerator_)];
        }
    }
    double structuralMutationProb = 0.02;
    if (probDist(randomGenerator_) < structuralMutationProb && mutatedIndividual.processExecutionAttemptSequence.size() > 1) {
         std::uniform_int_distribution<size_t> posDist(0, mutatedIndividual.processExecutionAttemptSequence.size());
         size_t pos = posDist(randomGenerator_);
         if (probDist(randomGenerator_) < 0.5 && mutatedIndividual.processExecutionAttemptSequence.size() > static_cast<size_t>(minSequenceLength_)) {
              if (pos < mutatedIndividual.processExecutionAttemptSequence.size()){
                   mutatedIndividual.processExecutionAttemptSequence.erase(mutatedIndividual.processExecutionAttemptSequence.begin() + pos);
              }
         } else {
              std::string randomProcess = availableProcessNames_[procIdxDist(randomGenerator_)];
               if (pos > mutatedIndividual.processExecutionAttemptSequence.size()) {
                    pos = mutatedIndividual.processExecutionAttemptSequence.size();
               }
              mutatedIndividual.processExecutionAttemptSequence.insert(mutatedIndividual.processExecutionAttemptSequence.begin() + pos, randomProcess);
         }
    }
    return mutatedIndividual;
}

// selectNextGeneration (inchangé)
void GeneticAlgorithm::selectNextGeneration() {
    std::vector<Individual> newPopulation;
    newPopulation.reserve(populationSize_);
    std::sort(currentPopulation_.begin(), currentPopulation_.end());
    for (int i = 0; i < eliteCount_ && i < static_cast<int>(currentPopulation_.size()); ++i) {
        newPopulation.push_back(currentPopulation_[i]);
    }
    std::vector<size_t> parentIndices = selectParentsViaRouletteWheel();
    if(parentIndices.empty() && populationSize_ > eliteCount_) {
         std::cerr << "Warning: No parents selected (empty indices), filling with elites or randoms." << std::endl;
         while(newPopulation.size() < static_cast<size_t>(populationSize_)) {
             if (!currentPopulation_.empty()) newPopulation.push_back(currentPopulation_[0]);
             else newPopulation.push_back(createRandomIndividual());
         }
         currentPopulation_ = newPopulation;
         return;
    }
     if (parentIndices.empty()) {
          currentPopulation_ = newPopulation;
          return;
     }

    std::uniform_int_distribution<size_t> parentDist(0, parentIndices.size() - 1);
    while (newPopulation.size() < static_cast<size_t>(populationSize_)) {
        size_t idx1_idx = parentDist(randomGenerator_);
        size_t idx2_idx = parentDist(randomGenerator_);
        if (idx1_idx >= parentIndices.size() || idx2_idx >= parentIndices.size()) { continue; }
        size_t index1 = parentIndices[idx1_idx];
        size_t index2 = parentIndices[idx2_idx];
         if (index1 >= currentPopulation_.size() || index2 >= currentPopulation_.size()) { continue; }

        int tries = 0;
        while (index1 == index2 && parentIndices.size() > 1 && tries < 10) {
            idx2_idx = parentDist(randomGenerator_);
             if (idx2_idx >= parentIndices.size()) continue;
            index2 = parentIndices[idx2_idx];
             if (index2 >= currentPopulation_.size()) continue;
            tries++;
        }
        const Individual& parent1 = currentPopulation_[index1];
        const Individual& parent2 = currentPopulation_[index2];
        std::pair<Individual, Individual> children = performCrossover(parent1, parent2);
        Individual child1 = applyMutation(children.first);
        Individual child2 = applyMutation(children.second);
        if (newPopulation.size() < static_cast<size_t>(populationSize_)) { newPopulation.push_back(child1); }
        if (newPopulation.size() < static_cast<size_t>(populationSize_)) { newPopulation.push_back(child2); }
    }
    currentPopulation_ = newPopulation;
}

// getBestIndividual (inchangé)
Individual GeneticAlgorithm::getBestIndividual() const {
    if (currentPopulation_.empty()) { return Individual(); }
    auto bestIt = std::max_element(currentPopulation_.begin(), currentPopulation_.end(),
                                   [](const Individual& a, const Individual& b) { return a.fitnessScore < b.fitnessScore; });
    return *bestIt;
}

// runEvolution (inchangé)
Individual GeneticAlgorithm::runEvolution(int numberOfGenerations) {
    std::cout << "Starting genetic algorithm evolution for " << numberOfGenerations << " generations..." << std::endl;
    initializePopulation();
    Individual bestOverallIndividual;
    // Initialiser bestOverallIndividual avec une fitness très basse pour la première comparaison
    bestOverallIndividual.fitnessScore = std::numeric_limits<double>::lowest();

    for (int gen = 0; gen < numberOfGenerations; ++gen) {
        evaluatePopulationFitness();
        Individual currentBest = getBestIndividual();
        // Mettre à jour le meilleur global si nécessaire
        if (gen == 0 || currentBest.fitnessScore > bestOverallIndividual.fitnessScore) {
             bestOverallIndividual = currentBest;
             std::cout << "Generation " << (gen + 1) << "/" << numberOfGenerations
                       << " - New Best Fitness: " << bestOverallIndividual.fitnessScore << std::endl;
        } else {
              std::cout << "Generation " << (gen + 1) << "/" << numberOfGenerations
                       << " - Best Fitness: " << currentBest.fitnessScore
                       << " (Overall Best: " << bestOverallIndividual.fitnessScore << ")" << std::endl;
        }
        selectNextGeneration();
    }
    evaluatePopulationFitness();
    Individual finalGenBest = getBestIndividual();
     if (finalGenBest.fitnessScore > bestOverallIndividual.fitnessScore) {
          bestOverallIndividual = finalGenBest;
     }
    std::cout << "Evolution finished. Final best fitness found: " << bestOverallIndividual.fitnessScore << std::endl;
    return bestOverallIndividual;
}
