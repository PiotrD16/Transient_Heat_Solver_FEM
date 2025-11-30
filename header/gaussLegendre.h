#pragma once
#include <vector>
#include <cmath>
#include <iostream>

using std::cout, std::vector;

// GaussLegendreSchema bedzie przechowywac wagi i punkty dla podanej liczby punktow

struct GaussLegendreSchema {
  private: vector < double > weights;
  vector < double > nodes;
  int elements;

  public: GaussLegendreSchema(const vector < double > & nodes, const vector < double > & weights):
      weights(weights), nodes(nodes), elements(static_cast < int > (nodes.size())) {};

  [[nodiscard]]const vector < double > getWeights() const {return weights;}
  [[nodiscard]]const vector < double > getNodes() const {return nodes;}
  [[nodiscard]]const int getElements() const {return elements;}
};

struct Integral {
private:
  GaussLegendreSchema schema1N;
  GaussLegendreSchema schema2N;
  GaussLegendreSchema schema3N;
  GaussLegendreSchema schema4N;
public:
  Integral():
      schema1N({0.0}, {2.0}),
      schema2N({-1.0 / std::sqrt(3.0),1.0 / std::sqrt(3.0)}, {1.0,1.0}),
      schema3N({-std::sqrt(3.0 / 5.0),0.0,std::sqrt(3.0 / 5.0)}, {5.0 / 9.0,8.0 / 9.0,5.0 / 9.0}),
      schema4N({
        -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)), // approx -0.861136
        -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)), // approx -0.339981
        std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),  // approx  0.339981
        std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0))   // approx  0.861136
      }, {
        // Corresponding Weights, in the exact same order
        (18.0 - std::sqrt(30.0)) / 36.0, // Weight for -0.861136 node
        (18.0 + std::sqrt(30.0)) / 36.0, // Weight for -0.339981 node
        (18.0 + std::sqrt(30.0)) / 36.0, // Weight for  0.339981 node
        (18.0 - std::sqrt(30.0)) / 36.0  // Weight for  0.861136 node
      }) {};
  ~Integral() = default;

  [[nodiscard]] const GaussLegendreSchema getSchema1N() const {return schema1N;}
  [[nodiscard]] const GaussLegendreSchema getSchema2N() const {return schema2N;}
  [[nodiscard]] const GaussLegendreSchema getSchema3N() const {return schema3N;}
  [[nodiscard]] const GaussLegendreSchema getSchema4N() const {return schema4N;}

  double calculateIntegralOneDim(double (*func)(double), int points) const;
  double calulateIntegralTwoDim(double (*func)(double ,double), int points) const;
};