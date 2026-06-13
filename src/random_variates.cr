# This library implements random variate generation for several
# common statistical distributions.
#
# Each distribution is implemented in its own class, and different
# parameterizations are created as instances of the class using the
# constructor to specify the parameterization.  All constructors use
# named parameters for clarity, so the order of parameters does not
# matter. All RV classes provide an optional argument *rng*, with which
# the user can specify a U(0,1) that is a subclass of *Random* to use as
# the core source of randomness. If *rng* is not specified, it defaults
# to *Random*. The *Random* class is extended below to implement a *next*
# method, making it (and subclasses) iterable.
#
# Once a random variate class has been instantiated, values can either be
# generated on demand using the *next* method or by using the instance as
# a generator in any iterable context.
#
module RandomVariates
  VERSION = "0.1.0"

  # Generate values uniformly distributed between *min* and *max*.
  #
  # *Arguments*::
  #   - *min* -> the lower bound for the range (default: 0).
  #   - *max* -> the upper bound for the range (default: 1).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Uniform
    include Iterator(Float64)

    @range : Float64

    getter :min, :max, :range

    def initialize(@min = 0.0, @max = 1.0, @rng : Random = Random.new)
      raise "Max must be greater than min." if max <= min
      @range = @max - @min
    end

    def next : Float64
      @min + @range * @rng.next
    end
  end

  # Triangular random variate generator with specified *min*, *mode*, and *max*.
  #
  # *Arguments*::
  #   - *min* -> the lower bound for the range.
  #   - *max* -> the upper bound for the range.
  #   - *mode* -> the highest likelihood value (*min* ≤ *mode* ≤ *max*).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Triangle
    include Iterator(Float64)

    getter :min, :max, :mode
    @range : Float64
    @crossover_p : Float64

    def initialize(@min = 0.0, @max = 1.0, @mode = 0.5, @rng : Random = Random.new)
      @range = @max - @min
      raise "Min must be less than Max." if @range <= 0
      raise "Mode must be between Min and Max." unless (@min..@max).includes? @mode

      @crossover_p = (@mode - @min) / @range
    end

    def next : Float64
      u = @rng.next
      u < @crossover_p ? @min + Math.sqrt(@range * (@mode - @min) * u) : @max - Math.sqrt(@range * (@max - @mode) * (1.0 - u))
    end
  end

  # Exponential random variate generator with specified *rate* or *mean*.
  # One and only one of *rate* or *mean* should be specified.
  #
  # *Arguments*::
  #   - *rate* -> the rate of occurrences per unit time (default: *nil*).
  #   - *mean* -> the expected value of the distribution (default: *nil*).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Exponential
    include Iterator(Float64)

    @mean : Float64 = -1.0
    @rate : Float64 = -1.0

    getter :rate, :mean

    def initialize(*, rate : Float64? = nil, mean : Float64? = nil, @rng : Random = Random.new)
      raise "Supply one and only one of mean or rate" unless rate.nil? ^ mean.nil?
      unless mean.nil?
        @mean = mean
        @rate = 1.0 / mean
      end
      unless rate.nil?
        @rate = rate
        @mean = 1.0 / rate
      end
      raise "Rate/mean must be strictly positive." if @rate <= 0
    end

    def next : Float64
      -@mean * Math.log(@rng.next)
    end
  end

  # Gaussian/normal random variate generator with specified
  # *mu* and +sigma deviation+.  Defaults to a standard normal.
  #
  # *Arguments*::
  #   - *mu* -> the expected value (default: 0).
  #   - *sigma* -> the standard deviation (default: 1).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Normal
    include Iterator(Float64)

    GAUSS_BOUND = 2.0 * Math.sqrt(2.0 / Math::E)

    getter :mu, :sigma

    def initialize(@mu = 0.0, @sigma = 1.0, @rng : Random = Random.new)
      raise "Sigma must be strictly positive." if @sigma <= 0.0
    end

    def next : Float64
      # Ratio of Uniforms
      loop do
        u = @rng.next
        next if u == 0.0
        v = GAUSS_BOUND * (@rng.next - 0.5)
        x = v / u
        x_sqr = x * x
        u_sqr = u * u
        if 6.0 * x_sqr <= 44.0 - 72.0 * u + 36.0 * u_sqr - 8.0 * u * u_sqr
          return @sigma * x + @mu
        elsif x_sqr * u >= 2.0 - 2.0 * u_sqr
          next
        elsif x_sqr <= -4.0 * Math.log(u)
          return @sigma * x + @mu
        end
      end
    end
  end

  # Alternate Gaussian/normal random variate generator with specified
  # *mu* and *sigma*.  Defaults to a standard normal.
  #
  # *Arguments*::
  #   - *mu* -> the expected value (default: 0).
  #   - *sigma* -> the standard deviation (default: 1).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class BoxMuller
    include Iterator(Float64)

    getter :mu, :sigma

    def initialize(@mu = 0.0, @sigma = 1.0, @rng : Random = Random.new)
      raise "Standard deviation must be positive." if @sigma <= 0
      @next_z = 0.0
      @need_new_pair = false
    end

    def next : Float64
      loop do
        @need_new_pair ^= true
        if @need_new_pair
          u = @rng.next
          v = @rng.next
          theta = 2.0 * Math::PI * u
          d = @sigma * Math.sqrt(-2.0 * Math.log(v))
          @next_z = @mu + d * Math.sin(theta)
          return @mu + d * Math.cos(theta)
        else
          return @next_z
        end
      end
    end
  end

  # Gamma generator based on Marsaglia and Tsang method Algorithm 4.33
  #
  # Produces gamma RVs with expected value *alpha* * *beta*.
  #
  # *Arguments*::
  #   - *alpha* -> the shape parameter (*alpha* > 0; default: 1).
  #   - *beta* -> the rate parameter (*beta* > 0; default: 1).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Gamma
    include Iterator(Float64)

    getter :alpha, :beta

    def initialize(@alpha = 1.0, @beta = 1.0, @rng : Random = Random.new)
      raise "Alpha and beta must be strictly positive." if alpha <= 0 || beta <= 0
      @std_normal = Normal.new(rng: @rng)
    end

    def next : Float64
      __gen__(@alpha, @beta)
    end

    private def __gen__(alpha, beta)
      if alpha > 1
        z = v = 0.0
        d = alpha - 1.0 / 3.0
        c = (1.0 / 3.0) / Math.sqrt(d)
        loop do
          loop do
            z = @std_normal.next
            v = 1.0 + c * z
            break if v > 0
          end
          z2 = z * z
          v = v * v * v
          u = @rng.next
          break if (u < 1.0 - 0.0331 * z2 * z2) ||
                   (Math.log(u) < (0.5 * z2 + d * (1.0 - v + Math.log(v))))
        end
        d * v * beta
      else
        result = __gen__(alpha + 1.0, beta)
        result * (@rng.next**(1.0 / alpha))
      end
    end
  end

  # Weibull generator based on Devroye
  #
  # *Arguments*::
  #   - *rate* -> the scale parameter (*rate* > 0; default: 1).
  #   - *k* -> the shape parameter (*k* > 0; default: 1).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Weibull
    include Iterator(Float64)

    getter :rate, :k

    def initialize(@rate = 1.0, @k = 1.0, @rng : Random = Random.new)
      raise "Rate and k must be positive." if @rate <= 0 || @k <= 0
      @power = 1.0 / @k
    end

    def next
      (-Math.log(@rng.next))**@power / @rate
    end
  end

  # Erlang generator - Weibull restricted to integer *k*
  #
  # *Arguments*::
  #   - *rate* -> the scale parameter (*rate* > 0; default: 1).
  #   - *k* -> the shape parameter (*k* > 0; default: 1).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Erlang < Weibull
    include Iterator(Float64)

    getter :rate, :k

    def initialize(@rate = 1.0, @k = 1, @rng : Random = Random.new)
      raise "K must be integer." unless k.integer?
      super(rate: @rate, k: @k, rng: @rng)
    end

    def next
      super.next
    end
  end

  # von Mises generator.
  #
  # This von Mises distribution generator is based on the VML algorithm by
  # L. Barabesis: "Generating von Mises variates by the Ratio-of-Uniforms Method"
  # Statistica Applicata Vol. 7, #4, 1995
  # http://sa-ijas.stat.unipd.it/sites/sa-ijas.stat.unipd.it/files/417-426.pdf
  #
  # *Arguments*::
  #   - *kappa* -> concentration coefficient (*kappa* ≥ 0).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class VonMises
    include Iterator(Float64)

    getter :kappa

    def initialize(@kappa, @rng : Random = Random.new)
      raise "kappa must be positive." if @kappa < 0
      @s = (@kappa > 1.3 ? 1.0 / Math.sqrt(@kappa) : Math::PI * Math.exp(-@kappa))
    end

    def next
      loop do
        r1 = @rng.next
        theta = @s * (2.0 * rng.next - 1.0) / r1
        next if theta.abs > Math::PI
        return theta if (
                          0.25 * @kappa * theta * theta < 1.0 - r1
                        ) || (
                          0.5 * @kappa * (Math.cos(theta) - 1.0) >= Math.log(r1)
                        )
      end
    end
  end

  # Poisson generator.
  #
  # *Arguments*::
  #   - *rate* -> expected number per unit time/distance (*rate* > 0; default: 1).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Poisson
    include Iterator(UInt32)

    @threshold : Float64

    getter :rate

    def initialize(@rate : Float64 = 1.0, @rng : Random = Random.new)
      raise "rate must be strictly positive." if @rate <= 0.0
      @threshold = Math.exp(-rate)
    end

    # Allow change of rate, which is much less computationally heavy
    # than instantiating new objects for every possible rate.
    def rate=(rate : Float64)
      raise "rate must be strictly positive." if @rate <= 0.0
      @threshold = Math.exp(-rate)
      @rate = rate
    end

    def next : UInt32
      count : UInt32 = 0
      product = 1.0
      loop do
        product *= @rng.next
        return count if product < @threshold
        count += 1
      end
    end
  end

  # Geometric generator.  Number of trials until first "success".
  #
  # *Arguments*::
  #   - *p* -> the probability of success (0 < *p* < 1; default: 0.5).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Geometric
    include Iterator(UInt32)

    @log_q : Float64

    getter :p

    def initialize(@p = 0.5, @rng : Random = Random.new)
      raise "Require 0 < p < 1." if @p <= 0 || @p >= 1
      @log_q = Math.log(1.0 - @p)
    end

    def next : UInt32
      (Math.log(1.0 - @rng.next) / @log_q).ceil.to_u
    end
  end

  # Binomial generator.  Number of "successes" in *n* independent trials.
  #
  # *Arguments*::
  #   - *n* -> the number of trials (*n* > 0, integer; default: 1).
  #   - *p* -> the probability of success (0 < *p* < 1; default: 0.5).
  #   - *rng* -> the (`Iterable`) source of U(0, 1)'s (default: `Random.new`)
  #
  class Binomial
    include Iterator(UInt32)

    @log_q : Float64

    getter :n, :p

    def initialize(@n : UInt32 = 1, @p = 0.5, @rng : Random = Random.new)
      raise "N must be a strictly positive integer." if @n <= 0
      raise "Require 0 < p < 1." if @p <= 0 || @p >= 1
      @complement = false
      if @p <= 0.5
        @log_q = Math.log(1.0 - @p)
      else
        @log_q = Math.log(@p)
        @complement = true
      end
    end

    def next : UInt32
      x : UInt32 = 0
      sum = 0.0
      loop do
        sum += Math.log(@rng.next) / (@n - x)
        return (@complement ? @n - x : x) if sum < @log_q
        x += 1
      end
    end
  end
end

module Random
  include Iterator(Float64)

  # Extends the base *Random* module to provide an *Iterator* based method for
  # generating U(0,1) values, which are the core for generating other distributions.
  def next
    next_float
  end
end
