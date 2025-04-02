#pragma once

#include <limits>
#include <string>
#include <vector>
#include <numeric> // Pour std::numeric_limits

// Représente un individu dans la population de l'algorithme génétique.
// Il contient une séquence de noms de processus à tenter et son score de fitness.
class Individual {
  public:
    // Constructeur par défaut (fitness très basse)
    Individual() : fitnessScore(std::numeric_limits<double>::lowest()) {}

    // Constructeur avec une séquence initiale
    explicit Individual(const std::vector<std::string>& sequence)
        : processExecutionAttemptSequence(sequence),
          fitnessScore(std::numeric_limits<double>::lowest()) {}

    // La séquence de noms de processus que le simulateur essaiera d'exécuter.
    // Ce n'est PAS un plan rigide, mais une liste de priorités/suggestions.
    std::vector<std::string> processExecutionAttemptSequence;

    // Le score de fitness calculé par le simulateur pour cette séquence.
    // Plus le score est élevé, meilleure est la solution.
    double fitnessScore;

    // Opérateur de comparaison pour permettre le tri des individus (utile pour l'élitisme).
    // Trie du meilleur (fitness élevé) au moins bon (fitness bas).
    bool operator<(const Individual& other) const {
        // Note: Tri décroissant pour que le meilleur soit au début après sort()
        return fitnessScore > other.fitnessScore;
    }
     // Opérateur de comparaison > pour priority_queue ou autres besoins
    bool operator>(const Individual& other) const {
        return fitnessScore < other.fitnessScore;
    }
};
