# Transit Delay Management

Tristan Ford

*The University of British Columbia*

## Costs

### Vehicles

This [study](https://www.columbia.edu/~ja3041/Electric%20Bus%20Analysis%20for%20NYC%20Transit%20by%20J%20Aber%20Columbia%20University%20-%20May%202016.pdf) from Columbia University, completed in 2016, identifies the lifetime cost of a transit vehicle to be $1.348M US which is nearly $1.9M CAD.  Therefore, using a conservative estimate of $2M CAD as the lifetime cost of purchasing a transit vehicle, we can calculate the ‘daily’ cost of adding a vehicle to the fleet.

Assuming a lifetime of 12 years (as is consistent with the Columbia study), by investing $2M CAD at a 5% rate of return, one could obtain $8,477 CAD biweekly until the investment is paid out.  This comes out to roughly $600 CAD per day, which we use as our cost for adding a new vehicle into the optimization.

### Operational Costs

The [BC Transit Service Plan](https://www.bctransit.com/wp-content/uploads/215/749/bct0.pdf) notes an estimated cost of operating conventional transit service as roughly $161 CAD per hour in 24/25.  Thus, we use a value of $160 CAD per hour in our model.  This is the cost of non-productive time spent between delivering service.  Furthermore, we translate deadhead distance to time by a fixed average speed of 30 km/hr and apply the same cost to penalize deadheading.

### Passenger Wait Time

The VTPI recommends a waiting time be valued at 100% of wages in this [report](https://www.vtpi.org/tca/tca0502.pdf).  In June 2024, the average hourly wage in BC was $36.63 as shown in this [report](https://www2.gov.bc.ca/assets/gov/data/statistics/people-population-community/income/earnings_and_employment_trends_data_tables.pdf) published by the government of BC.  Therefore, we use an hourly cost of $37 for passenger waiting time in our optimization.

### Return to Depot

TODO

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

Each trip delay distribution is given a mean value between 0-10% of its trip length, chosen randomly.  The standard deviation is set between 0-50% of its trip length, chosen randomly.  Both upper limits can be specified and we can set a random seed to guarantee reproducible results.

#### Ridership

We assume a uniform distribution between 0-120 passengers as the potential ridership for any given trip.  Each trip in an instance is assigned a value randomly between 0-120 and remains the same for each delay scenario.

The key difference is that ridership is not stochastic in the current model, though it could be in the future.

## Training/Testing

### Stochastic Delays

We demonstrate the value of considering stochastic delays by analyzing VSP instances with a range of trips.  We create 100 scenarios for each instance with different delays by sampling the distribution assigned to each trip.  All scenarios are optimized individually and we show the mean passenger delay for the optimized schedules in light blue.  We then calculate the mean passenger delay for each schedule with the delays from the other 99 scenarios and plot this in orange.  The plot on the bottom shows the similarity of the optimized schedules to the minimum cost solution for each instance.

<p align="center">
  <img src="./imgs/VSP-PD-200-100_1-scenario_comparison.png " />
</p>

We observe that as the number of trips increases, the model overfits to its training data and obtains very low passenger delays.  The performance on the test data degrades as a result.

Furthermore, the similarity with the minimum cost solution appears to settle to a ~40%.  This may indicate that some links that exist in the minimum cost solution are common across all solutions.  It would be interesting to see if this was the case and how to identify these links, which could serve as a basis for a heuristic method to improve solution times.

### Robustness to Delay

Using a series of instances with 30, 60, and 100 trips, we incrementally increase the variance of the delay distribution for each trip.  The variance increments are prescribed by setting the standard deviation equal to percentages of trip length and range from 10-100%.  For each variance increment, we optimize the minimum cost + delay model (represented in blue) over 100 scenarios.  We then plot the mean passenger delay and its standard deviation, considering a new ‘test’ set of 100 delay scenarios.

<p align="center">
  <img src="./imgs/VSP-PD-30-100-std_spectrum.png " />
</p>

<p align="center">
  <img src="./imgs/VSP-PD-60-100-std_spectrum.png " />
</p>

<p align="center">
  <img src="./imgs/VSP-PD-100-100-std_spectrum.png " />
</p>

Overall, we see the mean passenger delay and its standard deviation are significantly lower for the minimum delay + cost model than that of the minimum cost only.  Interestingly, the 60 trip test shows a lower passenger delay overall than the 30 and 100 trip tests.  This could be due to the ridership distribution and this may be evidence that the ridership profile is a critical factor in the optimization.

The similarities appear to hover around 60% and typically decrease as delays worsen.  Finding the common links between these solutions may yield heuristic methods to improve solution times as well as highlight 'manageial insights' regarding trip scheduling.

## Metrics

### Similarity

We calculate the similarity of two schedules using the `compareSchedules` function.  This returns the intersection over the union of the adjacency matrices for each schedule.  Schedules which share the same routing of vehicles between trips will have higher similarity.