require 'matrix'

module Stochastic
  class Fbm
    attr_accessor :hurst, :length, :delta

    def initialize(hurst, length, delta = 1.0)
      @hurst = hurst
      @length = length
      @delta = delta
    end

    def sample(hurst = @hurst, length = @length, delta = @delta)
      # Generate the covariance matrix
      cov_matrix = Matrix.build(length, length) do |i, j|
        delta * (0.5 * ((i+1)**(2*hurst) + (j+1)**(2*hurst) - (i-j).abs**(2*hurst)))
      end

      # Calculate the Cholesky decomposition of the covariance matrix
      cholesky_matrix = cholesky_factor(cov_matrix)

      # Generate a vector of Gaussian random variables
      rand_vector = Matrix.build(length, 1) { randn }

      # Calculate the fractional Brownian motion path
      fbm_path = cholesky_matrix * rand_vector

      return fbm_path.column(0).to_a
    end

    private
    def randn
      # Box-Muller method for generating normally distributed random variables
      theta = 2 * Math::PI * rand
      rho = Math.sqrt(-2 * Math.log(1 - rand))
      scale = rho * Math.cos(theta)

      return scale
    end

    def cholesky_factor(matrix)
      raise ArgumentError, "must provide symmetric matrix" unless matrix.symmetric?
      row_size = matrix.row_size
      l = Array.new(row_size) {Array.new(row_size, 0)}
      (0 ... row_size).each do |k|
        (0 ... row_size).each do |i|
          if i == k
            sum = (0 .. k-1).inject(0.0) {|sum, j| sum + l[k][j] ** 2}
            val = Math.sqrt(matrix[k,k] - sum)
            l[k][k] = val
          elsif i > k
            sum = (0 .. k-1).inject(0.0) {|sum, j| sum + l[i][j] * l[k][j]}
            val = (matrix[k,i] - sum) / l[k][k]
            l[i][k] = val
          end
        end
      end
      Matrix[*l]
    end
  end
end
