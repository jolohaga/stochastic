## Description

Stochastic - A collection of stochastic processes.

### `Stochastic::Fbm`

The `Stochastic::Fbm` function takes three arguments:

1. The Hurst parameter (`hurst`)
2. The length of the time series (`length`)
3. A `delta` parameter that scales the increments.  Optional

It generates a fractional Brownian motion using the covariance matrix and Cholesky decomposition method.

The covariance matrix is created based on the Hurst parameter, and the Cholesky decomposition is performed on this matrix.

A vector of Gaussian random variables is generated using the randn function, which uses the Box-Muller method.

The fractional Brownian motion path is computed as the product of the Cholesky matrix and the vector of random variables.

The path is then returned as an array.

#### Usage

To use this script, set your Hurst parameter and time length, call the `Stochastic::Fbm#sample` function, and use or print the generated path as needed.

    hurst_parameter = 0.75
    time_length = 1000
    fbm_path = Stochastic::Fbm.new(hurst_parameter, time_length)
    
    # Print or use the fbm_path
    fbm_path.sample.each_with_index do |val, index|
      puts "t#{index}: #{val}"
    end
