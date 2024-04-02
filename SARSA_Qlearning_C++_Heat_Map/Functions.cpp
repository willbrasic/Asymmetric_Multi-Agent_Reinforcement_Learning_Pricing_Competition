#include "Functions.h"


//////////////////////////////
// Function Definitions
//////////////////////////////


// Function to populate each agent's action space
std::vector<double> populate_A(
    int m,
    double A_min,
    double step_size 
)
{   
    // Initialize vector to store action space for each agent
    std::vector<double> action_space(m, 0.0);

    // Populate action space
    for (int j{}; j < m; ++j)
    {
        action_space[j] = A_min + (j * step_size);
    }

    // Return the action space
    return action_space;
}

// Cartesian product of action space
std::vector<std::vector<double>> generate_combinations(
    const std::vector<double>& A,
    int n,
    int m
)
{
    // All possible combinations of actions for agents
    std::vector<std::vector<double>> As;
    for (int i{ 0 }; i < n; ++i)
    {
        if (i == 0)
        {
            for (int j = 0; j < m; ++j)
            {
                As.push_back({ A[j] });
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
                    combination.push_back(A[j]);
                    newAs.push_back(combination);
                }
            }
            As = newAs;
        }
    }

    return As;
}

// Function to determine quantity for each state for each agent
std::vector<std::vector<double>> q_fn(
    const std::vector<std::vector<double>>& price,
    const std::vector<int>& a,
    int n,
    int a_0,
    double mu,
    int S_cardinality
)
{
    // First obtain the demoninator of demand function
    std::vector<std::vector<double>> quantity_numerator(S_cardinality, std::vector<double>(n, 0.0));
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            quantity_numerator[j][i] = exp((a[i] - price[j][i]) / mu);
        }
    }
    // Finish computing demoninator of quantity
    std::vector<double> quantity_denominator(S_cardinality, 0.0);
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        quantity_denominator[j] = std::accumulate(quantity_numerator[j].begin(),quantity_numerator[j].end(),exp(a_0 / mu));
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


// Function to determine total market revenue for each state for each agent
std::vector<double> rvn_fn(
    const std::vector<std::vector<double>>& price,
    const std::vector<std::vector<double>>& quantity,
    int S_cardinality,
    int n
)
{
    // Initialize vector to store revenue for each state for each agent
    std::vector<std::vector<double>> revenue_each_agent(S_cardinality, std::vector<double>(n, 0.0));
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            revenue_each_agent[j][i] = price[j][i] * quantity[j][i];
        }
    }

    // Initialize vector to store total market revenue for each state
    std::vector<double> revenue(S_cardinality, 0.0);
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        revenue[j] = std::accumulate(revenue_each_agent[j].begin(), revenue_each_agent[j].end(), 0.0);
    }

    return revenue;
}

// Function to determine profit
std::vector<std::vector<double>> r_fn(
        const std::vector<std::vector<double>>& price,
        const std::vector<std::vector<double>>& quantity,
        const std::vector<int>& marginal_cost,
        int S_cardinality,
        int n
)
{
    // Initialize vector to store profits
    std::vector<std::vector<double>> profit(S_cardinality, std::vector<double>(n, 0.0));
    for (int j{ 0 }; j < S_cardinality; ++j)
    {
        for (int i{ 0 }; i < n; ++i)
        {
            profit[j][i] = (price[j][i] - marginal_cost[i]) * quantity[j][i];
        }
    }

    return profit;
}

// Function to determine consumer surplus for each state
std::vector<double> cs_fn(
    const std::vector<std::vector<double>>& price,
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
            cs_1[j][i] = exp((a[i] - price[j][i]) / mu);
        }
    }

    // Sum state consumer surplus across agents
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


// Function to determine average profits for each agent when setting each of the m possible prices
std::vector<std::vector<double>> Qi0_fn(
    const std::vector<std::vector<double>>& price,
    const std::vector<std::vector<double>>& profit,
    const std::vector<double>& A,
    int n,
    int m,
    double delta
)
{
    // Initialize matrix for storing average profits for each agent for each possible price
    std::vector<std::vector<double>> Qi0(m, std::vector<double>(n, 0.0));

    // Loop over agents
    for (int i { 0 }; i < n; ++i)
    {
        // Loop over prices
        for (int j { 0 }; j < m; ++j)
        {
            // Initialize vector to store indices of states where agent i sets the j'th price
            std::vector<int> indices;

            // Loop over the possible states
            for (int k{ 0 }; k < price.size(); ++k)
            {
                // If in state k agent i sets price j
                if (price[k][i] == A[j])
                {
                    // Add this state index to indices vector
                    indices.push_back(k);
                }
            }

            // Calculate profit summation for agent i when setting price j
            double sum{ 0.0 };
            for (int index : indices) {
                sum += profit.at(index)[i];
            }
            // Find number of times agent i set price j
            int count = indices.size();

            // Calculate mean profit for agent i when setting price j (avoiding division by 0)
            Qi0[j][i] = (count > 0) ? (sum / count) / (1.0 - delta) : 0.0;
        }
    }
    return Qi0;
} 

// Function to populate Q-matrix for each agent
std::vector<std::vector<std::vector<double>>> Q_fn(
    const std::vector<std::vector<double>>& Qi0,
    int n,
    int m,
    int S_cardinality
)
{
    // Initialize Q-matrix for storage
    std::vector<std::vector<std::vector<double>>> Q(S_cardinality, std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0)));

    // Initialize each agent's matrix with average profits when setting each of the m possible prices
    for (int i{ 0 }; i < n; ++i)
    {   
        // Loop over number of states
        for (int s{ 0 }; s < S_cardinality; ++s)
        {   
            // Loop of number of prices
            for (int p{ 0 }; p < m; ++p)
            {   
                // Initialize Q-matrix in episode e for agent i when setting price p in state s as the 
                // average profit when agent i sets price p
                Q[s][p][i] = Qi0[p][i];
            }
        }
    }

    return Q;
}

// Function to generate all possible values of alpha
std::vector<double> alpha_fn(
        int alpha_values,
        double alpha_min,
        double alpha_step_size
)
{
    // Create vector to store values of alpha
    std::vector<double> alpha_vector(alpha_values, 0.0);

    // Populate alpha_vector with alpha_values equally spaced points
    for (int alpha{ 0 }; alpha < alpha_values; ++alpha)
    {
        alpha_vector[alpha] = alpha_min + (alpha * alpha_step_size);
    }

    return alpha_vector;
}

// Function to generate all possible values of beta
std::vector<double> beta_fn(
        int beta_values,
        double beta_min,
        double beta_step_size
)
{
    // Create vector to store values of beta
    std::vector<double> beta_vector(beta_values, 0.0);

    // Populate the beta_vector with beta_values equally spaced points
    for (int beta{ 0 }; beta < beta_values; ++beta)
    {
        beta_vector[beta] = beta_min + (beta * beta_step_size);
    }

    return beta_vector;
}



