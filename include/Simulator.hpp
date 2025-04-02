#pragma once

#include "Config.hpp" // Contient les définitions de Process, Stock, etc.
#include "Process.hpp"
#include <map>
#include <string>
#include <vector>
#include <utility> // Pour std::pair
#include <queue>   // Pour std::priority_queue

// Simule l'exécution des processus en fonction d'une séquence guide,
// en gérant le temps, les ressources, le parallélisme et les contraintes.
// Calcule le score de fitness d'une séquence donnée.
class Simulator {
  public: // Rendre public pour que ProcessManager puisse utiliser RunningProcessInfo

    // Structure interne pour suivre un processus en cours d'exécution.
    // Déplacée ici pour être publique.
    struct RunningProcessInfo {
        const Process* processDetails; // Pointeur vers la définition du processus
        int completionTime;            // Cycle auquel le processus se termine

        // Opérateur de comparaison pour la file de priorité (min-heap).
        // Permet de récupérer le processus qui se termine le plus tôt.
        bool operator>(const RunningProcessInfo& other) const {
            return completionTime > other.completionTime;
        }
    };

  public: // Interface publique du simulateur
    // Constructeur
    Simulator(const Config& config, int timeLimit);

    // Fonction CLÉ : Évalue une séquence et retourne son fitness.
    // Gère la simulation temporelle, les ressources, le parallélisme.
    // Remplit également les logs d'exécution si la simulation réussit.
    // @param sequenceToAttempt: La séquence de noms de processus à essayer (le "chromosome").
    // @param executionLogs: Vecteur de sortie pour stocker les paires {temps_démarrage, nom_processus}.
    // @return: Le score de fitness calculé pour cette séquence.
    double calculateFitnessAndLogs(
        const std::vector<std::string>& sequenceToAttempt,
        std::vector<std::pair<int, std::string>>& executionLogs);

  private:
    // Référence à la configuration globale (processus, stocks initiaux, objectifs).
    const Config& systemConfig_;
    // Limite de temps pour la simulation (le "delay" fourni en argument).
    int simulationTimeLimit_;

    // Helper : Vérifie si un processus peut démarrer avec les stocks actuels.
    // @param process: Le processus à vérifier.
    // @param currentStocks: L'état actuel des stocks.
    // @return: true si le processus peut démarrer, false sinon.
    bool canProcessStart(const Process* process, const std::map<std::string, int>& currentStocks) const;

    // Helper : Trouve un processus dans la configuration par son nom.
    // @param name: Le nom du processus à rechercher.
    // @return: Un pointeur constant vers le processus trouvé, ou nullptr si non trouvé.
    const Process* findProcessByName(const std::string& name) const;
};
