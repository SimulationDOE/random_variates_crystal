# random_variates

This library implements random variate generation for several
common statistical distributions.

Each distribution is implemented in its own class, and different
parameterizations are created as instances of the class using the
constructor to specify the parameterization.  All constructors use
named parameters for clarity, so the order of parameters does not
matter. All *RandomVariate* classes provide an optional argument *rng*,
with which the user can specify a U(0,1) that is a subclass of *Random*
to use as the core source of randomness. If *rng* is not specified, it
defaults to *Random*. The *Random* class is extended to implement
a *next* method, making it (and alternative random number generators
written as subclasses) iterable.

Once a random variate class has been instantiated, values can either be
generated on demand using the *next* method or by using the instance as
a generator in any iterable context.


## Installation

1. Add the dependency to your `shard.yml`:

   ```yaml
   dependencies:
     random_variates:
       github: SimulationDOE/random_variates_crystal
   ```

2. Run `shards install`

## Usage

```crystal
require "random_variates"

tri = RandomVariates::Triangle.new(min: -5.0, max: 5.0, mode: 3.0)
pois = RandomVariates::Poisson.new(rate: 2.0)
exp = RandomVariates::Exponential.new(mean: 42)  # specify either mean or rate

# generate one-by-one with next
5.times { puts "#{tri.next},#{pois.next},#{exp.next}" }

puts # to get a separator blank line

# generate an array with 3 standard normals
puts RandomVariates::Normal.new.first(3).to_a
```

Generator classes are available for Uniform, Triangle, Exponential, Normal,
Gamma, Poisson, Geometric, and Binomial distributions. Details are available
in the docs or the source code.

## Contributors

- [SimulationDOE](https://github.com/SimulationDOE/random_variates_crystal) - creator and maintainer
