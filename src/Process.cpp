#include "Process.hpp"
#include "Stock.hpp"

Process::Process(const std::string& name, int cycleDuration)
    : name_(name), cycleDuration_(cycleDuration) {}

void Process::addRequirement(const std::string& stockName, int quantity) {
    requirements_[stockName] = quantity;
}

void Process::addResult(const std::string& stockName, int quantity) {
    results_[stockName] = quantity;
}

const std::string& Process::getName() const {
    return name_;
}

int Process::getCycleDuration() const {
    return cycleDuration_;
}

bool Process::canExecute(const std::unordered_map<std::string, std::shared_ptr<Stock>>& stocks) const {
    for (const auto& [stockName, requiredQuantity] : requirements_) {
        auto it = stocks.find(stockName);
        if (it == stocks.end() || it->second->getQuantity() < requiredQuantity) {
            return false;
        }
    }
    return true;
}

bool Process::execute(std::unordered_map<std::string, std::shared_ptr<Stock>>& stocks) {
    // First check if we have enough resources
    if (!canExecute(stocks)) {
        return false;
    }
    
    // Consume required stocks
    for (const auto& [stockName, requiredQuantity] : requirements_) {
        stocks[stockName]->remove(requiredQuantity);
    }
    
    // Produce result stocks
    for (const auto& [stockName, producedQuantity] : results_) {
        auto it = stocks.find(stockName);
        if (it != stocks.end()) {
            it->second->add(producedQuantity);
        } else {
            // Create the stock if it doesn't exist
            stocks[stockName] = std::make_shared<Stock>(stockName, producedQuantity);
        }
    }
    
    return true;
}

const std::unordered_map<std::string, int>& Process::getRequirements() const {
    return requirements_;
}

const std::unordered_map<std::string, int>& Process::getResults() const {
    return results_;
}