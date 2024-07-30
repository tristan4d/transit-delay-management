# Transit Delay Management

Tristan Ford

*The University of British Columbia*

## Costs

### Vehicles

[Aber, J. (2016)][1] identifies the lifetime cost of a transit vehicle to be $1.348M US which is nearly $1.9M CAD.  Therefore, using a conservative estimate of $2M CAD as the lifetime cost of purchasing a transit vehicle, we can calculate the ‘daily’ cost of adding a vehicle to the fleet.

Assuming a lifetime of 12 years (as is consistent with the Columbia study), by investing $2M CAD at a 5% rate of return, one could obtain $8,477 CAD biweekly until the investment is paid out.  This comes out to roughly $600 CAD per day, which we use as our cost for adding a new vehicle into the optimization.

### Operational Costs

The BC Transit [Service Plan (2023)][2] notes an estimated cost of operating conventional transit service as roughly $161 CAD per hour in 24/25.  Thus, we use a value of $160 CAD per hour in our model.  This is the cost of any time spent for a vehicle between leaving and returning to the depot.  Furthermore, we translate deadhead distance to time by a fixed average speed of 30 km/hr and apply the same cost to penalize deadheading.

### Passenger Wait Time

[Litman, T. (2016)][3] recommends a waiting time be valued at 100% of wages.  In June 2024, the average hourly wage in BC was $36.63 as reported by [BC Stats (2024)][4].  Therefore, we use an hourly cost of $37 for passenger waiting time in our optimization.

### Pull Out/In Trips

The distance from the depot to/from the start/end of a vehicle schedule is costed at the operational cost level by transforming the distance to time using an average speed of 30 km/hr.  These are called pull out and pull in trips, respectively.

When two trips are separated by 3 hours or more, it is preferable to return to the depot in between as time spent at the depot can be valued at $0.  This also helps in reducing unnecessary increases in vehicle numbers in the model as a singular trip that is far in time from others in a given vehicle schedule can be incorporated without causing undue cost and making use of the same vehicle.

## Data

### Historical

#### Trips

To build vehicle scheduling problem instances, we follow the network flow model wherein each trip is represented by a node in a network.  Two trips are said to be compatible if both may be feasibly operated by the same vehicle, meaning the start and end locations and times are such that one vehicle could operate one trip first and the other afterwards.  Arcs are created between compatible trips in the direction of vehicle flows.  This creates a directed acyclic graph (ignoring the depot) as our instance of the VSP.

To create this graph, we rely upon General Transit Feed Specification (GTFS) data, which is a widely used format for documenting and sharing transit information.  We are able to obtain trip information (including start and stop locations and times) from these files which allows us to create the input graph to the VSP.

#### Delays

TODO: when received from BC Transit

#### Ridership

TODO: when received from BC Transit

### Generated

#### Delays

Each trip will have a distribution of delays that we must understand to properly model the minimization of network delays.  We assume a normal distribution for each trip centered around zero, as the trip length ought to be increased if the average trip delay is much greater than zero.  We can then select a standard deviation as a percentage of each trip’s length to assign a ‘realistic’ delay distribution to each individual trip.  Then, when creating scenarios, we sample each trip’s distribution to obtain a snapshot of the network for one planning horizon.

Each trip delay distribution is given a mean value between -10 to +10% of its trip length, chosen randomly.  The standard deviation is set between 10-50% of its trip length, chosen randomly.  Both upper limits can be specified and we can set a random seed to guarantee reproducible results.

#### Ridership

We assume a uniform distribution between 0-120 passengers as the potential ridership for any given trip.  Each trip in an instance is assigned a value randomly between 0-120 and remains the same for each delay scenario.

The key difference is that ridership is not stochastic in the current model, though it could be in the future.

#### Depot

We are able to specify a depot location `(lat, lon)` and if none is specified, we select the mean starting location of all trips as the depot location.  This somewhat mimics the preference for agencies to position their depots closest to where most of their trips operate.

## Why Use Stochastic Delays?

### Stochastic Delays

We demonstrate the value of considering stochastic delays by analyzing VSP instances with a range of trips.  We create 10 different instances for each trip increment from 25-100.  For each trip increment, we optimize three models.  The first is the minimum cost solution to the instance, disregarding delay information.  The sceond is a minimum cost + delay model optimized over the sample mean delay of each trip only.  The third is a minimum cost + delay model optimized over stochastic delays sampled from each trip distribution.  The results are calculated against an unseen sample of delays as to replicate optimizing over historical delay information and performance on real time disruptions.

![overfit_unseen](imgs/VSP-PD-100-10_10-overfit_unseen.png)

We observe that as the number of trips increases, both the minimum cost and mean trip delay models overfits to its training data and the delays that occur when considering the unseen sample grow.  The performance on the stochastically optimized model, however, maintains low delays on both the mean trip delay and stochastic delay scenarios.

We also observe that the costs for the minimum cost and mean trip delay models are quite similar, indicating that the minimum cost solution is sufficient if one were only concerned with the mean trip delays in a system with reasonable on-time-performance.  The stochastically optimized model exhibits lower passenger delays on the unseen data, as expected, as well lower cost.  This highlights the models overall superiority when considering network disruptions.

Furthermore, the similarity to the minimum cost solution decreases in both delay considering models.  This makes sense as when we add trips, the number of feasible schedules grows and the solution which minimizes some delay is likely to be different than that of the minimum cost solution.  Lastly, we note that the utilization, as compared to the minimum cost solution, is nearly identical and this suggests these models are not adding much buffer time in the vehicle schedules.

### Robustness to Delay

Using an instance with 100 trips, we incrementally increase the variance of the delay distribution for each trip.  The variance increments are prescribed by setting the standard deviation equal to percentages of trip length and range from 10-100%.  For each variance increment, we optimize three models.  The first is the minimum cost solution to the instance, disregarding delay information.  The sceond is a minimum cost + delay model optimized over the sample mean delay of each trip only.  The third is a minimum cost + delay model optimized over stochastic delays sampled from each trip distribution.  The results are calculated against an unseen sample of delays as to replicate optimizing over historical delay information and performance on real time disruptions.


![delay_spectrum_unseen](imgs/VSP-PD-100-10_10-delay_spectrum_unseen.png)

We note that the minimum cost solutions and the model optimized over sample mean delays only perform very similarly.  This makes sense as the mean delays are small and therefore the model does not need to adjust the schedules much to accommodate these minor disruptions.  Interestingly, as the variance increases, we see these two models have a high similarity of ~80%.  When the variance is low, the sampled scenarios will be closer to the mean of each distribution and therefore the sample mean is likely to be positive.  However, when the variance is high and because we truncate the distributions at a maximum value of 1x trip length, the sample means are likely to be negative which means all trips arrive early.  This is exactly the same as finding the minimum cost solution which may explain why we see this high similarity at large variances.

The model optimized over the stochastic delay sample performs much better, when considering delays on any given day.  Even when delay distributions have high variance, the average passenger delay per trip and the cost of service are kept low.  We see the similarity with the minimum cost solution degrades and the utilization dips slightly.  This suggests that the stochstically optimized model is selecting higher cost links with more buffer time as this is less expensive than the large passenger delays that are occuring in the other solutions.

An interesting artifact from this plot is in the performance of the two delay-considering models on the sample mean instance.  We see that the stochastically optimized model has lower passenger delays even though the other model was trained specifically on this data.  The overall cost for the model trained on sample means is lower than the stochastically optimized model, though, which indicates the former is selecting cheaper links at the expense of causing more passenger delays.  This tradeoff results in lower costs when considering the sample mean delays only, but we observe how this impacts individual scenarios through the stochastic cost.

### How Many Scenarios?

We demonstrate that a relatively small number of scenarios are required to obtain stable results in the stochastically optimized model.  We create 10 different instances with 100 trips for each number of scenario increment from 1-10.  For each trip increment, we optimize three models.  The first is the minimum cost solution to the instance, disregarding delay information.  The sceond is a minimum cost + delay model optimized over the sample mean delay of each trip only.  The third is a minimum cost + delay model optimized over stochastic delays sampled from each trip distribution.  The results are calculated against an unseen sample of delays as to replicate optimizing over historical delay information and performance on real time disruptions.

![scenarios_unseen](imgs/VSP-PD-100-10_1-10-scenarios_unseen.png)

We observe that by ~10 different scenarios, the stochastically optimized model has settled to a steady state in terms of passenger delays, cost, and similarity to the minimum cost solution.  For instances with 30 trips and 50 trips, we saw settling occur at ~3 and ~5 scenarios, respectively.  This suggests that the solution is not changing drastically as we add more scenarios beyond ~10% of the trip total.  As a result, we set the number of scenarios for any instance to 10% of its number of trips.

## Metrics

### Utilization

TODO

### Passenger Delay

TODO

### Cost

TODO

### Similarity

We calculate the similarity of two schedules using the `compareSchedules` function.  This returns the intersection over the union of the adjacency matrices for each schedule.  Schedules which share the same routing of vehicles between trips will have higher similarity.

## Comparison to Other Solutions

### Run Time Analysis

Run time analysis (RTA) consists of analyzing historical trip run times and adjusting the planned travel time to align with current conditions.  Common practice is to set the planned travel time to allow a given percentile of historical trips to have been made on time - typically a value of 85% is used.  In this section, we simulate performing a RTA on each instance by increasing travel times for all trips to accommodate the 85th percentile of delay scenarios.  We then re-optimize this newly created instance as a minimum cost flow model.

First, we consider how these models perform when the variance of trip delays increases.  We allow trip delay distributions to have a maximum variance of 20-100% of their trip length.  We expect the stochastically optimized model to perform the best, in terms of cost, as a large variance will mean many trips are also arriving early as well as very late.  Therefore, the RTA model will likely overcompensate by adding significant run time to each trip which will increase cost.

![rta_std_comparison](imgs/VSP-PD-100-10_10-rta_std_comparison.png)

Indeed, we see that the stochasically optimized model outperforms the other two in terms of cost.  The RTA model does outperform the others with respect to passenger delays, which makes sense as the travel times for each trip are longer.  In this analysis, a key assumption is that the delay distributions would remain the same, regardless of trip length, meaning we simply shift the mean of each trip delay distribution by the amount of travel time added to that trip.  This is a large assumption.

Next, we consider how these models perform when the mean of trip delays increases.  Similarly, we allow the means to be larger, but we expect the RTA to outperform the others in terms of cost as delays become more extreme.  This is because when delays become extreme across the entire system, it may not be feasible to arrange the trips in a way that avoids large passenger delays.

![rta_mean_comparison](imgs/VSP-PD-100-10_10-rta_mean_comparison.png)

Interestingly, we see the stochastically optimized model still has the lowest cost.

**TODO**

* excess run time costs may be what is missing
* looking at average trip passenger delays may be watering down delays - should look at inf norm?

## Miscellaneous

### Stochastic Delay Vignette

Consider the mathematical model for our problem,

$$
\begin{gathered}
\min&\sum_{i,j\in V}c_{ij}x_{ij}+\sum_{k\in\mathcal{S},i\in T}r_is_i^k& \\
\text{subject to:}&x_{ij}\in\mathbb{Z}_+&\forall i,j,\in V \\
&s_i^k\geq0&\forall i\in T,k\in\mathcal{S} \\
&s_i^k\geq\sum_{j\in T}x_{ji}(s_j^k+l_j^k-b_{ji})&\forall i\in T,k\in\mathcal{S} \\
&\sum_{j\in V}x_{ij}=\sum_{j\in V}x_{ji}&\forall i\in V \\
&\sum_{j\in V}x_{ij}=1&\forall i\in T,
\end{gathered}
$$

where $s_i^k$ is defined as the amount of time that trip $i$ departs after its scheduled departure time in scenario $k$.  Note, this cannot be negative as trips are not allowed to depart early.  $l_i^k$ is the primary delay for trip $i$ in scenario $k$.

Suppose the mean primary delay for all trips is less than or equal to 0, that is

$$
\frac{1}{|\mathcal{S}|}\sum_{k\in\mathcal{S}}l_i^k\leq0\quad\forall i\in T.
$$

Should we restrict our model to the scenario including only mean primary delays, then we observe that $s_i^k=0$ is optimal and the model collapses to the minimum cost solution.  However, should we sample the distribution of each trip to build multiple delay scenarios, we will inevitably encouter scenarios with positive primary delays (assuming a normal distribution).  As we have established, the minimum cost solution can perform significantly worse than one which considers delays in terms of delay propagation and cost.  Thus, the naive approach of optimizing over the mean primary delay for all trips is not recommended for developing schedules that are robust to system disruptions.

## References

[1]: https://www.columbia.edu/~ja3041/Electric%20Bus%20Analysis%20for%20NYC%20Transit%20by%20J%20Aber%20Columbia%20University%20-%20May%202016.pdf "Aber, J. (2016). Electric bus analysis for new york city transit. Columbia University, 3."
[2]: https://www.bctransit.com/wp-content/uploads/215/749/bct0.pdf "BC Transit (2023). 2023/24 - 2025/26 service plan."
[3]: https://www.vtpi.org/tca/tca0502.pdf "Litman, T. (2016). Transportation cost and benefit analysis II - travel time costs 5.2 travel time and speed."
[4]: https://www2.gov.bc.ca/assets/gov/data/statistics/people-population-community/income/earnings_and_employment_trends_data_tables.pdf "BC Stats (2024). Earnings & employment trends - june 2024."