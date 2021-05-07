#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <mav_trajectory_generation/polynomial_optimization_linear.h>
#include <mav_trajectory_generation/polynomial_optimization_nonlinear.h>


using namespace Eigen;

// see https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix
template<typename M>
M load_csv (const std::string & path) {
  std::ifstream indata;
  indata.open(path);
  std::string line;
  std::vector<double> values;
  uint rows = 0;
  while (std::getline(indata, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ',')) {
      values.push_back(std::stod(cell));
    }
    ++rows;
  }
  return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}

int main(int argc, char **argv)
{
  std::string inputFile;
  std::string outputFile;
  double v_max;
  double a_max;

  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("input,i", po::value<std::string>(&inputFile)->required(), "path to list of waypoints")
    ("output,o", po::value<std::string>(&outputFile)->required(), "path to output files")
    ("v_max", po::value<double>(&v_max)->default_value(1.0), "maximum velocity [m/s]")
    ("a_max", po::value<double>(&a_max)->default_value(1.0), "maximum velocity [m/s^2]")
  ;

  try
  {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 0;
    }
  }
  catch(po::error& e)
  {
    std::cerr << e.what() << std::endl << std::endl;
    std::cerr << desc << std::endl;
    return 1;
  }

  if (!boost::filesystem::exists(inputFile)) {
    std::cerr << "No such file: " << inputFile << std::endl;
    return 1;
  }

  mav_trajectory_generation::Vertex::Vector vertices;
  const int dimension = 3;
  const int derivative_to_optimize = mav_trajectory_generation::derivative_order::JERK;

  // create vertices with their constraints
  auto input = load_csv<Matrix<double, Dynamic, 3> >(inputFile);

  if (input.rows() < 2) {
    std::cerr << "Not enough datapoints given!" << std::endl;
    return 1;
  }

  if (input.rows() == 2) {
    mav_trajectory_generation::Vertex v1(dimension);
    v1.makeStartOrEnd(input.row(0), derivative_to_optimize);
    vertices.push_back(v1);

    mav_trajectory_generation::Vertex v2(dimension);
    auto middle = (input.row(1) - input.row(0)) * 0.5 + input.row(0);
    v2.addConstraint(mav_trajectory_generation::derivative_order::POSITION, middle);
    vertices.push_back(v2);

    mav_trajectory_generation::Vertex v3(dimension);
    v3.makeStartOrEnd(input.row(1), derivative_to_optimize);
    vertices.push_back(v3);
  } else {
    // at least 3 points given
    for (int row = 0; row < input.rows(); ++row) {
      // std::cout << input.row(row) << std::endl;
      mav_trajectory_generation::Vertex vertex(dimension);
      if (row == 0 || row == input.rows() - 1) {
        vertex.makeStartOrEnd(input.row(row), derivative_to_optimize);
      } else {
        vertex.addConstraint(mav_trajectory_generation::derivative_order::POSITION, input.row(row));
      }
      vertices.push_back(vertex);
    }
  }

  // compute segment times
  std::vector<double> segment_times;
  const double magic_fabian_constant = 6.5; // A tuning parameter.
  segment_times = estimateSegmentTimes(vertices, v_max, a_max, magic_fabian_constant);

  // solve
  const int N = 8;
  mav_trajectory_generation::Segment::Vector segments;

#if SOLVE_LINEAR
  mav_trajectory_generation::PolynomialOptimization<N> opt(dimension);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  opt.solveLinear();

  // Obtain the polynomial segments.
  opt.getSegments(&segments);
#else
  mav_trajectory_generation::NonlinearOptimizationParameters parameters;
  parameters.max_iterations = 1000;
  parameters.f_rel = 0.05;
  parameters.x_rel = 0.1;
  parameters.time_penalty = 500.0;
  parameters.initial_stepsize_rel = 0.1;
  parameters.inequality_constraint_tolerance = 0.1;

  mav_trajectory_generation::PolynomialOptimizationNonLinear<N> opt(dimension, parameters, false);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  opt.addMaximumMagnitudeConstraint(mav_trajectory_generation::derivative_order::VELOCITY, v_max);
  opt.addMaximumMagnitudeConstraint(mav_trajectory_generation::derivative_order::ACCELERATION, a_max);
  opt.optimize();

  // Obtain the polynomial segments.
  opt.getPolynomialOptimizationRef().getSegments(&segments);
#endif

  std::ofstream output(outputFile);

  output << "Duration,x^0,x^1,x^2,x^3,x^4,x^5,x^6,x^7,y^0,y^1,y^2,y^3,y^4,y^5,y^6,y^7,z^0,z^1,z^2,z^3,z^4,z^5,z^6,z^7,yaw^0,yaw^1,yaw^2,yaw^3,yaw^4,yaw^5,yaw^6,yaw^7" << std::endl;
  Eigen::IOFormat csv_fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", ",", "", "", "", "");

  for (const auto& segment : segments) {
    output << segment.getTime() << ",";
    for (const auto& polynomial : segment.getPolynomialsRef()) {
      Eigen::VectorXd coefficients = polynomial.getCoefficients();
      output << coefficients.format(csv_fmt) << ",";
    }
    output << "0,0,0,0,0,0,0,0" << std::endl;
  }

  return 0;
}
