# KRPSim - Process Chain Optimizer

KRPSim is a process chain optimization tool that finds optimal execution schedules for interconnected processes.

## Project Structure

The project consists of two main programs:

1. **krpsim**: The optimizer program that finds the best solution for executing processes
2. **krpsim_verif**: A verification program that checks if a solution is valid

## Building the Project

To build the project, simply run:

```
make
```

This will create the `krpsim` and `krpsim_verif` executables in the `bin` directory.

## Usage

### krpsim

```
./bin/krpsim <file> <delay>
```

- `<file>`: The configuration file with stocks and processes
- `<delay>`: The maximum time (in seconds) the program will run to find a solution

### krpsim_verif

```
./bin/krpsim_verif <file> <result_to_test>
```

- `<file>`: The configuration file with stocks and processes
- `<result_to_test>`: The file containing the solution to verify

## Configuration File Format

The configuration file format is as follows:

- Lines starting with `#` are comments
- Stocks are defined as `<stock_name>:<quantity>`
- Processes are defined as `<name>:(<need>:<qty>[;<need>:<qty>[...]]):(<result>:<qty>[;<result>:<qty>[...]]):duration`
- Optimization targets are defined as `optimize:(<stock_name>|time[;<stock_name>|time[...]])`

Example:
```
# stock name:quantity
euro:10

# process name:(needs):(results):delay
equipment_purchase:(euro:8):(equipment:1):10
product_creation:(equipment:1):(product:1):30
delivery:(product:1):(happy_client:1):20

# optimize
optimize:(time;happy_client)
```

## Solution File Format

The solution file format is simple:

```
<cycle>:<process_name>
```

For example:
```
0:equipment_purchase
10:product_creation
40:delivery
```

## How It Works

The optimizer uses a combination of greedy algorithms and local search to find a good solution. The approach is:

1. Parse the configuration file to get stocks, processes, and optimization targets
2. Generate an initial greedy solution
3. Improve the solution using local search
4. Return the best solution found within the time limit

The verification program:

1. Parses the configuration file and the solution file
2. Simulates the execution of the solution step by step
3. Checks that no stock goes negative
4. Reports whether the solution is valid and the final stock quantities

## Modern C++ Features Used

This project uses several modern C++ features:

- Smart pointers (`std::shared_ptr`, `std::unique_ptr`)
- Move semantics
- Standard library containers (`std::vector`, `std::unordered_map`)
- `std::optional` for nullable return values
- Range-based for loops
- Auto type deduction
- Lambda expressions
- Regular expressions
- Chrono library for time measurement

## Contributors

This project is a solution to the "krpsim" algorithmic project.