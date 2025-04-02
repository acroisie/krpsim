#pragma once

#include "Config.hpp"
#include "Individual.hpp"
#include "Simulator.hpp"
#include <vector>
#include <string>
#include <random> // Pour la génération aléatoire

// Implémente l'algorithme génétique pour trouver une séquence d'actions
// quasi-optimale pour le problème de planification KrpSim.
class GeneticAlgorithm {
  public:
    // Constructeur
    // @param config: La configuration chargée (processus, stocks initiaux, objectifs).
    // @param simulator: Le simulateur utilisé pour évaluer les individus.
    // @param populationSize: Le nombre d'individus dans chaque génération.
    // @param mutationRate: La probabilité qu'une mutation survienne sur un gène (processus).
    // @param crossoverRate: La probabilité qu'un croisement survienne entre deux parents.
    // @param eliteCount: Le nombre des meilleurs individus à conserver directement (élitisme).
    // @param minSequenceLength: Longueur minimale pour les séquences initiales.
    // @param maxSequenceLength: Longueur maximale pour les séquences initiales.
    GeneticAlgorithm(const Config& config,
                       Simulator& simulator,
                       int populationSize = 100,         // Taille de population par défaut
                       double mutationRate = 0.05,       // Taux de mutation par défaut
                       double crossoverRate = 0.7,       // Taux de croisement par défaut
                       int eliteCount = 2,               // Nombre d'élites par défaut
                       int minSequenceLength = 50,       // Longueur min initiale
                       int maxSequenceLength = 150);     // Longueur max initiale

    // Exécute l'algorithme génétique pour un nombre donné de générations.
    // @param numberOfGenerations: Le nombre de cycles d'évolution à effectuer.
    // @return: Le meilleur individu trouvé après toutes les générations.
    Individual runEvolution(int numberOfGenerations);

    // Retourne le meilleur individu de la population actuelle.
    // NOTE: Peut être const si on recalcule/retrouve le meilleur sans trier sur place.
    //       Pour l'instant, on la laisse const et on utilise std::max_element.
    Individual getBestIndividual() const;

  private:
    // Références et paramètres de l'AG
    const Config& systemConfig_;          // Configuration du problème
    Simulator& planningSimulator_;        // Simulateur pour l'évaluation
    int populationSize_;
    double mutationRate_;
    double crossoverRate_;
    int eliteCount_;
    int minSequenceLength_;
    int maxSequenceLength_;

    // État de l'AG
    std::vector<Individual> currentPopulation_; // La population actuelle d'individus
    std::vector<std::string> availableProcessNames_; // Cache des noms de processus valides

    // Générateur de nombres aléatoires pour l'AG
    // Doit être mutable car son état change à chaque génération de nombre.
    std::mt19937 randomGenerator_;

    // --- Fonctions internes du cycle de l'AG ---

    // Génère la population initiale avec des séquences aléatoires.
    void initializePopulation(); // Modifie currentPopulation_ et randomGenerator_ -> Non-const

    // Évalue le fitness de chaque individu dans la population actuelle
    // en utilisant le simulateur.
    void evaluatePopulationFitness(); // Modifie les fitness dans currentPopulation_ -> Non-const

    // Crée la population de la génération suivante via sélection, croisement et mutation.
    void selectNextGeneration(); // Modifie currentPopulation_ et randomGenerator_ -> Non-const

    // --- Opérateurs Génétiques ---

    // Sélectionne des individus parents en utilisant la méthode de la roue de la fortune.
    // Retourne les indices des parents sélectionnés dans currentPopulation_.
    // *** DOIT être non-const car utilise randomGenerator_ ***
    std::vector<size_t> selectParentsViaRouletteWheel();

    // Effectue un croisement (ici, double point) entre deux parents pour créer deux enfants.
    // *** DOIT être non-const car utilise randomGenerator_ ***
    std::pair<Individual, Individual> performCrossover(const Individual& parent1, const Individual& parent2);

    // Applique des mutations à un individu (ponctuelle et/ou structurelle).
    // *** DOIT être non-const car utilise randomGenerator_ ***
    Individual applyMutation(const Individual& individual);

    // Crée un nouvel individu avec une séquence aléatoire (pour l'initialisation).
    // *** DOIT être non-const car utilise randomGenerator_ ***
    Individual createRandomIndividual();
};
