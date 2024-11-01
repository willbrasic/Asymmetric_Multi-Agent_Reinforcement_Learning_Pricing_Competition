
#include "Functions.h"


//////////////////////////////
// Function Definitions
//////////////////////////////


// Function to populate each firms action space
std::vector<double> populate_action_space_fn(
        int m,
        double A_min,
        double A_step_size
)
{
    // Initialize vector to store action space for each firm
    std::vector<double> action_space(m, 0.0);

    // Populate action space
    for (int j{ 0 }; j < m; ++j)
    {
        action_space[j] = A_min + (j * A_step_size);
    }

    // Return the action space
    return action_space;
}

// Cartesian product of action space
std::vector<std::vector<double>> action_space_cartesian_product_fn(
        const std::vector<double>& action_space,
        int n,
        int m
)
{
    // All possible combinations of actions for firms
    std::vector<std::vector<double>> As;
    for (int i{ 0 }; i < n; ++i)
    {
        if (i == 0)
        {
            for (int j = 0; j < m; ++j)
            {
                As.push_back({ action_space[j] });
            }
        }
        else
        {
            std::vector<std::vector<double>> newAs;
            for (int j{ 0 }; j < m; ++j)
            {
                for (const auto& prev : As)
                {
                    std::vector<double> combination = prev;
                    combination.push_back(action_space[j]);
                    newAs.push_back(combination);
                }
            }
            As = newAs;
        }
    }

    return As;
}

// Function to determine quantity
std::vector<std::vector<double>> quantity_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<int>& a,
        int n,
        int a_0,
        double mu,
        int S_cardinality
)
{
    // First obtain the denominator of demand function
    std::vector<std::vector<double>> quantity_numerator(S_cardinality, std::vector<double>(n, 0.0));
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            quantity_numerator[j][i] = exp((a[i] - prices[j][i]) / mu);
        }
    }
    
    // Finish computing denominator of quantity
    std::vector<double> quantity_denominator(S_cardinality, 0.0);
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        quantity_denominator[j] = std::accumulate(
                quantity_numerator[j].begin(),
                quantity_numerator[j].end(),
                exp(a_0 / mu));
    }

    // Compute quantity
    std::vector<std::vector<double>> quantity(S_cardinality, std::vector<double>(n, 0.0));
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            quantity[j][i] = quantity_numerator[j][i] / quantity_denominator[j];
        }
    }

    return quantity;
}


// Function to determine total market revenue for each state
std::vector<double> revenue_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<std::vector<double>>& quantity,
        int S_cardinality,
        int n
)
{
    // Initialize vector to store revenue for each state for each firm
    std::vector<std::vector<double>> revenue_each_firm(S_cardinality, std::vector<double>(n, 0.0));
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            revenue_each_firm[j][i] = prices[j][i] * quantity[j][i];
        }
    }

    // Initialize vector to store total market revenue for each state
    std::vector<double> revenue(S_cardinality, 0.0);
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        revenue[j] = std::accumulate(revenue_each_firm[j].begin(), revenue_each_firm[j].end(), 0.0);
    }

    return revenue;
}

// Function to determine consumer surplus for each state
std::vector<double> cs_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<int>& a,
        double mu,
        int a_0,
        int n,
        int S_cardinality
)
{
    // Initialize vector to get first part of consumer surplus
    std::vector<std::vector<double>> cs_1(S_cardinality, std::vector<double>(n, 0.0));
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            cs_1[j][i] = exp((a[i] - prices[j][i]) / mu);
        }
    }

    // Sum state consumer surplus across firms
    std::vector<double> cs_1_sum(S_cardinality, 0.0);
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        cs_1_sum[j] = std::accumulate(cs_1[j].begin(), cs_1[j].end(), exp(a_0 / mu));
    }

    // Final consumer surplus
    std::vector<double> consumer_surplus(S_cardinality, 0.0);
    for (int j { 0 }; j < S_cardinality; ++j)
    {
        consumer_surplus[j] = mu * log(cs_1_sum[j]);
    }

    return consumer_surplus;
}

// Function to determine profits
std::vector<std::vector<double>> profits_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<std::vector<double>>& quantity,
        const std::vector<int>& marginal_cost,
        int S_cardinality,
        int n
)
{
    // Initialize vector to store profits
    std::vector<std::vector<double>> profits(S_cardinality, std::vector<double>(n, 0.0));
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            profits[j][i] = (prices[j][i] - marginal_cost[i]) * quantity[j][i];
        }
    }

    return profits;
}

// Function to determine average profits for both firms from setting each possible price
std::vector<std::vector<double>> avg_profits_for_states_fn(
        const std::vector<std::vector<double>>& prices,
        const std::vector<std::vector<double>>& profits,
        const std::vector<double>& action_space,
        int n,
        int m,
        double delta
)
{
    // Initialize matrix for storing average profitss for each firm for each possible price
    std::vector<std::vector<double>> avg_profits_for_states_matrix(m, std::vector<double>(n, 0.0));

    // Loop over firms
    for (int i { 0 }; i < n; ++i)
    {
        // Loop over prices
        for (int j { 0 }; j < m; ++j)
        {
            // Initialize vector to store indices of states where firm i sets the j'th price
            std::vector<int> indices;

            // Loop over the possible states
            for (int k{ 0 }; k < prices.size(); ++k)
            {
                // If in state k firm i sets prices j
                if (prices[k][i] == action_space[j])
                {
                    // Add this state index to indices vector
                    indices.push_back(k);
                }
            }

            // Calculate profits summation for firm i when setting price j
            double sum{ 0.0 };
            for (const int& index : indices)
            {
                sum += profits[index][i];
            }
            
            // Find number of times firm i set price j
            int count = indices.size();

            // Calculate mean profits for firm i when setting price j (avoiding division by 0)
            avg_profits_for_states_matrix[j][i] = (count > 0) ? (sum / count) / (1.0 - delta) : 0.0;
        }
    }
    return avg_profits_for_states_matrix;
}

// Function to populate Q-matrix
std::vector<std::vector<double>> initial_Q_matrix_fn(
        const std::vector<std::vector<double>>& avg_profits_for_states_matrix,
        int n,
        int m,
        int S_cardinality
)
{
    // Initialize Q_matrix
    std::vector<std::vector<double>> initial_Q_matrix(S_cardinality, std::vector<double>(m, 0.0));

    // Populate Q-matrix
    for (int i{ 0 }; i < n; ++i)
    {
        // Loop over number of states
        for (int s{ 0 }; s < S_cardinality; ++s)
        {
            // Loop of number of actions
            for (int j{ 0 }; j < m; ++j)
            {
                // Initialize Q-matrix in state s when taking action j as the average profits when firm i sets price p
                initial_Q_matrix[s][j] = avg_profits_for_states_matrix[j][i];
            }
        }
    }
    return initial_Q_matrix;
}

// Function to get logit competitive equilibrium outcome
std::vector<std::vector<double>> logit_competitive_outcome_fn(
        const std::vector<int>& a,
        const std::vector<int>& mc,
        double mu,
        int n,
        int a_0,
        int iteration_limit
)
{
    // Initial prices for each firm
    std::vector<double> p_0(n, 2.5);

    // Store distances between old and new prices
    std::vector<double> norm_vector(n, 0.0);

    // Initialize prices for convergence loop to p_0
    std::vector<double> competitive_prices = p_0;

    // Dampening parameter for competitive price convergence
    double eta{ 0.5 };

    // Tolerance for convergence
    double tol{ 1e-8 };

    // Average and maximum distances between old and new prices
    double avg_norm{ 0.0 }, max_norm{ 0.0 };

    // Convergence iteration counter for competitive case
    int competitive_convergence{ 0 };

    // Print up to eight decimal places and turn of scientific notation
    std::cout << std::setprecision(8) << std::fixed;

    // Do-while loop to get competitive outcome
    do
    {
        // Increment iteration tracker
        ++competitive_convergence;

        // Set p_0 as the old prices
        p_0 = competitive_prices;

        // Denominator for competitive quantity
        std::vector<double> exp_quantities(n, 0.0);
        double quantity_denominator = exp(a_0 / mu);
        for (int i{ 0 }; i < n; ++i)
        {
            double price = competitive_prices[i];
            exp_quantities[i] = exp((a[i] - price) / mu);
            quantity_denominator += exp_quantities[i];
        }

        // Competitive quantity
        std::vector<double> competitive_quantities(n, 0.0);
        for (int i{ 0 }; i < n; ++i)
        {
            double quantity = exp_quantities[i] / quantity_denominator;
            competitive_quantities[i] = quantity;
        }

        // Competitive price
        for (int i{ 0 }; i < n; ++i)
        {
            competitive_prices[i] = (eta * p_0[i]) + ((1 - eta) * (mc[i] + (mu / (1 - competitive_quantities[i]) )));
        }

        // Compute the distance between each price vector element-wise
        std::vector<double> price_norm(n);
        for (size_t i{ 0 }; i < n; i++)
        {
            price_norm[i] = std::abs(competitive_prices[i] - p_0[i]);
        }

        // Obtain the maximum distance
        max_norm = *max_element(price_norm.begin(), price_norm.end());

        // Compute the average distance
        avg_norm = std::accumulate(price_norm.begin(), price_norm.end(), 0.0) / n;

        // Create vector of the max and average distances
        norm_vector.emplace_back(max_norm);
        norm_vector.emplace_back(avg_norm);

        // If at *iteration_limit*, probably not convering so quit the loop
        if (competitive_convergence == iteration_limit)
        {
            std::cout << "Probably not converging so quitting fixed point iteration." << std::endl;
            max_norm = avg_norm = tol;
        }

    } while (max_norm > tol && avg_norm > tol);

    // Print up to four decimal places and turn of scientific notation
    std::cout << std::setprecision(4) << std::fixed;

    // Denominator for competitive quantity
    std::vector<double> exp_quantities(n, 0.0);
    double quantity_denominator = exp(a_0 / mu);
    for (int i{ 0 }; i < n; ++i)
    {
        double price = competitive_prices[i];
        exp_quantities[i] = exp((a[i] - price) / mu);
        quantity_denominator += exp_quantities[i];
    }

    // Converged competitive quantity
    std::vector<double> competitive_quantities(n, 0.0);
    for (int i{ 0 }; i < n; ++i)
    {
        double quantity = exp_quantities[i] / quantity_denominator;
        competitive_quantities[i] = quantity;
    }

    // Converged competitive revenues
    std::vector<double> competitive_revenues(n, 0.0);
    for (int i{ 0 }; i < n; ++i)
    {
        competitive_revenues[i] = competitive_prices[i] * competitive_quantities[i];
    }

    // Converged competitive profits
    std::vector<double> competitive_profits(n, 0.0);
    for (int i{ 0 }; i < n; ++i)
    {
        competitive_profits[i] = (competitive_prices[i] - mc[i]) * competitive_quantities[i];
    }

    // Log of summation term in consumer surplus formula
    double exp_cs_sum = exp( a_0 / mu );
    for (int i{ 0 }; i < n; ++i)
    {
        exp_cs_sum += exp((a[i] - competitive_prices[i]) / mu);
    }
    double log_exp_cs_sum = log(exp_cs_sum);

    // Consumer surplus at time step t in episode e
    double cs = (mu * log_exp_cs_sum) + a_0 / mu;
    std::vector<double> cs_vector(1, 0.0);
    cs_vector.emplace_back(cs);

    // Add competitive results to vector
    std::vector<std::vector<double>> competitive_outcome { competitive_profits, competitive_quantities, competitive_prices,
                                                           competitive_revenues, cs_vector  };

    return competitive_outcome;

}

// Function to print main results after all episodes
void print_results(
        const std::vector<double>& avg_profits,
        const std::vector<double>& avg_prices,
        const std::vector<double>& avg_quantities,
        const std::vector<double>& avg_revenues,
        const std::vector<int>& a,
        const std::vector<int>& mc,
        double avg_cs,
        double average_seconds,
        double mu,
        int average_time_steps,
        int n,
        int a_0,
        int iteration_limit
)
{
    // Competitive Outcome
    std::vector<std::vector<double>> competitive_outcome = logit_competitive_outcome_fn(a, mc, mu, n, a_0, iteration_limit);
    std::vector<double> competitive_profits = competitive_outcome[0];
    std::vector<double> competitive_quantities = competitive_outcome[1];
    std::vector<double> competitive_prices = competitive_outcome[2];
    std::vector<double> competitive_revenues = competitive_outcome[3];
    double competitive_cs = competitive_outcome[4][1];

    std::cout << std::fixed << std::setprecision(4);

    std::cout << "\n****************************************************\n";
    std::cout << "***************** OVERALL RESULTS ******************\n";
    std::cout << "****************************************************\n";
    std::cout << "\nAverage number of time steps until convergence: " << average_time_steps << "\n";
    std::cout << "\nAverage seconds until convergence: " << average_seconds << "\n";

    std::cout << "                       -----------------------------------------------------\n";
    std::cout << "             Tot/Avg";
    for (int i{ 1 }; i <= n; ++i)
    {
        std::cout << std::setw(16) << i;
    }
    std::cout << "\n----------------------------------------------------------------------------\n";

    // Total average profit among all agents
    double total_profits = std::accumulate(avg_profits.begin(), avg_profits.end(), 0.0);
    std::cout << "\nProfit        " << total_profits;

    // Individual agent average profit
    for (const double& profit : avg_profits)
        std::cout << std::setw(17) << profit;

    // Total average demand along all agents
    double total_quantities = std::accumulate(avg_quantities.begin(), avg_quantities.end(), 0.0);
    std::cout << "\nDemand        " << total_quantities;

    // Individual agent average demand
    for (const double& quantity : avg_quantities)
        std::cout << std::setw(17) << quantity;

    // Average price over all agents
    double avg_price = std::accumulate(avg_prices.begin(), avg_prices.end(), 0.0) / n;
    std::cout << "\nPrice         " << avg_price;

    // Individual agent average price
    for (const double& price : avg_prices)
        std::cout << std::setw(17) << price;

    // Total average revenue among all agents
    double total_revenues = std::accumulate(avg_revenues.begin(), avg_revenues.end(), 0.0);
    std::cout << "\nRevenue       " << total_revenues;
    // Individual agent average revenue
    for (const double& revenue : avg_revenues)
        std::cout << std::setw(17) << revenue;

    // Average consumer surplus
    std::cout << "\nCS            " << avg_cs;

    std::cout << std::fixed << std::setprecision(2);

    std::cout << "\n\n\n****************************************************\n";
    std::cout << "**** PERCENTAGE CHANGE FROM COMPETITIVE OUTCOME ****\n";
    std::cout << "****************************************************\n";
    std::cout << "\n                                 Firms                 \n";
    std::cout << "                       -----------------------------------------------------\n";
    std::cout << "             Tot/Avg";
    for (int i{ 1 }; i <= n; ++i)
    {
        std::cout << std::setw(16) << i;
    }
    std::cout << "\n----------------------------------------------------------------------------\n";

    // Total profit percentage change from competitive outcome
    double total_competitive_profits = std::accumulate(competitive_profits.begin(), competitive_profits.end(), 0.0);
    double total_profits_percentage_change = 100 * ( (total_profits - total_competitive_profits) / total_competitive_profits );
    std::cout << "\nProfit        " << total_profits_percentage_change << "%";

    // Individual agent profit percentage change from competitive outcome
    for (int i{ 0 }; i < n; ++i)
        std::cout << std::setw(17) << 100 * ( (avg_profits[i] - competitive_profits[i]) / competitive_profits[i] ) << "%";

    // Total demand percentage change from competitive outcome
    double total_competitive_quantities = std::accumulate(competitive_quantities.begin(), competitive_quantities.end(), 0.0);
    double total_quantities_percentage_change = 100 * ( (total_quantities - total_competitive_quantities) / total_competitive_quantities );
    std::cout << "\nDemand       " << total_quantities_percentage_change << "%";

    // Individual agent demand percentage change from competitive outcome
    for (int i{ 0 }; i < n; ++i)
        std::cout << std::setw(17) << 100 * ( (avg_quantities[i] - competitive_quantities[i]) / competitive_quantities[i] ) << "%";

    // Average price percentage change from competitive outcome
    double avg_competitive_price = std::accumulate(competitive_prices.begin(), competitive_prices.end(), 0.0) / n;
    double avg_price_percentage_change = 100 * ( (avg_price - avg_competitive_price) / avg_competitive_price );
    std::cout << "\nPrice         " << avg_price_percentage_change << "%";

    // Individual agent price percentage change from competitive outcome
    for (int i{ 0 }; i < n; ++i)
        std::cout << std::setw(17) << 100 * ( (avg_prices[i] - competitive_prices[i]) / competitive_prices[i] ) << "%";

    // Total revenues percentage change from competitive outcome
    double total_competitive_revenues = std::accumulate(competitive_revenues.begin(), competitive_revenues.end(), 0.0);
    double total_revenues_percentage_change = 100 * ( (total_revenues - total_competitive_revenues) / total_competitive_revenues );
    std::cout << "\nRevenue        " << total_revenues_percentage_change << "%";
    // Individual agent revenue percentage change from competitive outcome
    for (int i{ 0 }; i < n; ++i)
        std::cout << std::setw(17) << 100 * ( (avg_revenues[i] - competitive_revenues[i]) / competitive_revenues[i] ) << "%";

    // Consumer surplus percentage change from competitive outcome
    double cs_percentage_change = 100 * ( (avg_cs - competitive_cs) / competitive_cs );
    std::cout << "\nCS            " << cs_percentage_change;
}


