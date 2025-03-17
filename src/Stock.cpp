#include "Stock.hpp"

Stock::Stock(const std::string& name, int quantity)
    : name_(name), quantity_(quantity) {}

const std::string& Stock::getName() const {
    return name_;
}

int Stock::getQuantity() const {
    return quantity_;
}

void Stock::add(int amount) {
    quantity_ += amount;
}

bool Stock::remove(int amount) {
    if (quantity_ >= amount) {
        quantity_ -= amount;
        return true;
    }
    return false;
}

void Stock::setQuantity(int quantity) {
    quantity_ = quantity;
}