#pragma once

#include "Config.hpp"
#include "Individual.hpp"
#include "Simulator.hpp"
#include "GeneticAlgorithm.hpp" // Inclure GeneticAlgorithm
#include <map>
#include <string>
#include <vector>
#include <memory> // Pour std::unique_ptr si besoin

// Gère l'orchestration globale : chargement, exécution de l'AG,
// et affichage de la meilleure solution trouvée.
class ProcessManager {
  public:
    // Constructeur: initialise le simulateur et l'algorithme génétique.
    // @param config: La configuration du problème.
    // @param delayLimit: La limite de temps pour la simulation.
    ProcessManager(const Config& config, int delayLimit);

    // Lance le processus complet : exécute l'AG et affiche les résultats.
    void run();

  private:
    const Config& config_; // Configuration (passée par référence)
    Simulator simulator_;  // Le simulateur (possédé par ProcessManager)
    GeneticAlgorithm geneticAlgorithm_; // L'algorithme génétique (possédé par ProcessManager)
    int delayLimit_; // Stocker la limite de temps

    // Génère la sortie formatée attendue par krpsim_verif
    // à partir de la meilleure solution trouvée par l'AG.
    // @param bestSolution: L'individu représentant la meilleure séquence trouvée.
    void generateOutput(const Individual& bestSolution);
};
