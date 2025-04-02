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

// evaluatePopulationFitness (inchangé)
void GeneticAlgorithm::evaluatePopulationFitness() {
    for (auto& individual : currentPopulation_) {
        std::vector<std::pair<int, std::string>> logs;
        individual.fitnessScore = planningSimulator_.calculateFitnessAndLogs(
            individual.processExecutionAttemptSequence, logs);
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
         if (indiv.fitnessScore < minFitness) minFitness = indiv.fitnessScore;
         if (indiv.fitnessScore > maxFitness) maxFitness = indiv.fitnessScore;
     }

     // Si tous les fitness sont identiques (et potentiellement très bas), la sélection par roue est inutile/impossible.
     if (minFitness == maxFitness) {
          std::cout << "Warning: All individuals have the same fitness (" << minFitness << "). Selecting parents randomly." << std::endl;
          selectedIndices.reserve(populationSize_);
          std::uniform_int_distribution<size_t> randDist(0, currentPopulation_.size() - 1);
          for(int i = 0; i < populationSize_; ++i) {
              selectedIndices.push_back(randDist(randomGenerator_));
          }
          return selectedIndices;
     }


     // *** CORRECTION ICI: S'assurer que le shift ne crée pas de problèmes avec des nombres énormes ***
     // Utiliser une approche qui évite l'addition directe si minFitness est trop extrême.
     // On peut normaliser les scores entre 0 et (max - min).
     double range = maxFitness - minFitness;
     if (range <= 0 || !std::isfinite(range)) { // Si range est invalide ou nul (devrait être couvert par le check précédent)
         range = 1.0; // Éviter division par zéro
     }

     std::vector<double> normalizedFitnessValues;
     normalizedFitnessValues.reserve(currentPopulation_.size());
     double totalNormalizedFitness = 0.0;

     for (const auto& indiv : currentPopulation_) {
         // Normaliser entre 0 et 1 (approximativement, ajouter epsilon pour éviter 0 strict)
         double normalized = (indiv.fitnessScore - minFitness) / range + 0.001;
         normalizedFitnessValues.push_back(normalized);
         totalNormalizedFitness += normalized;
     }


     if (totalNormalizedFitness <= 0.0 || !std::isfinite(totalNormalizedFitness)) {
         // Si même après normalisation ça échoue, sélection aléatoire
         std::cerr << "Warning: Total normalized fitness is non-positive or invalid (" << totalNormalizedFitness << "). Selecting parents randomly." << std::endl;
          selectedIndices.reserve(populationSize_);
          std::uniform_int_distribution<size_t> randDist(0, currentPopulation_.empty() ? 0 : currentPopulation_.size() - 1);
          for(int i = 0; i < populationSize_; ++i) {
               if (currentPopulation_.empty()) break;
              selectedIndices.push_back(randDist(randomGenerator_));
          }
          return selectedIndices;
     }

     // Distribution pour choisir un point sur la roue normalisée
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


// performCrossover (inchangé par rapport à v2)
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

// applyMutation (inchangé par rapport à v2)
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

// selectNextGeneration (inchangé par rapport à v2)
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
     // Vérifier si parentIndices est vide *avant* de créer la distribution
     if (parentIndices.empty()) {
          // Gérer le cas où la population est remplie uniquement par élitisme
          currentPopulation_ = newPopulation;
          return;
     }

    std::uniform_int_distribution<size_t> parentDist(0, parentIndices.size() - 1);
    while (newPopulation.size() < static_cast<size_t>(populationSize_)) {
        size_t idx1_idx = parentDist(randomGenerator_);
        size_t idx2_idx = parentDist(randomGenerator_);
        // Vérifier les limites au cas où parentIndices serait plus petit que prévu
        if (idx1_idx >= parentIndices.size() || idx2_idx >= parentIndices.size()) {
             std::cerr << "Error: Invalid parent index selected." << std::endl;
             continue; // Sauter cette itération
        }
        size_t index1 = parentIndices[idx1_idx];
        size_t index2 = parentIndices[idx2_idx];
        // Vérifier les limites de currentPopulation_
         if (index1 >= currentPopulation_.size() || index2 >= currentPopulation_.size()) {
              std::cerr << "Error: Parent index out of bounds for current population." << std::endl;
              continue; // Sauter cette itération
         }

        int tries = 0;
        while (index1 == index2 && parentIndices.size() > 1 && tries < 10) {
            idx2_idx = parentDist(randomGenerator_);
             if (idx2_idx >= parentIndices.size()) continue; // Vérif supplémentaire
            index2 = parentIndices[idx2_idx];
             if (index2 >= currentPopulation_.size()) continue; // Vérif supplémentaire
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

// getBestIndividual (inchangé par rapport à v2)
Individual GeneticAlgorithm::getBestIndividual() const {
    if (currentPopulation_.empty()) { return Individual(); }
    auto bestIt = std::max_element(currentPopulation_.begin(), currentPopulation_.end(),
                                   [](const Individual& a, const Individual& b) { return a.fitnessScore < b.fitnessScore; });
    return *bestIt;
}

// runEvolution (inchangé par rapport à v2)
Individual GeneticAlgorithm::runEvolution(int numberOfGenerations) {
    std::cout << "Starting genetic algorithm evolution for " << numberOfGenerations << " generations..." << std::endl;
    initializePopulation();
    Individual bestOverallIndividual;
    for (int gen = 0; gen < numberOfGenerations; ++gen) {
        evaluatePopulationFitness();
        Individual currentBest = getBestIndividual();
        if (gen == 0 || bestOverallIndividual.fitnessScore < -1e17 /* init value */ || currentBest.fitnessScore > bestOverallIndividual.fitnessScore) {
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
     if (bestOverallIndividual.fitnessScore < -1e17 || finalGenBest.fitnessScore > bestOverallIndividual.fitnessScore) {
          bestOverallIndividual = finalGenBest;
     }
    std::cout << "Evolution finished. Final best fitness found: " << bestOverallIndividual.fitnessScore << std::endl;
    return bestOverallIndividual;
}
