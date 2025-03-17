#pragma once

#include <string>

/**
 * @class Stock
 * @brief Represents a resource/stock in the system
 */
class Stock {
public:
    /**
     * @brief Constructor
     * @param name The name of the stock
     * @param quantity The initial quantity
     */
    Stock(const std::string& name, int quantity);
    
    /**
     * @brief Get the name of the stock
     * @return The stock name
     */
    const std::string& getName() const;
    
    /**
     * @brief Get the current quantity
     * @return The quantity
     */
    int getQuantity() const;
    
    /**
     * @brief Add to the current quantity
     * @param amount The amount to add
     */
    void add(int amount);
    
    /**
     * @brief Remove from the current quantity
     * @param amount The amount to remove
     * @return True if successful, false if insufficient quantity
     */
    bool remove(int amount);
    
    /**
     * @brief Set the quantity directly
     * @param quantity The new quantity
     */
    void setQuantity(int quantity);
    
private:
    std::string name_;
    int quantity_;
};