![GitHub license](https://img.shields.io/github/license/torresflo/Simple-Math-Libray.svg)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](http://makeapullrequest.com)
![GitHub contributors](https://img.shields.io/github/contributors/torresflo/Simple-Math-Library.svg)
![GitHub issues](https://img.shields.io/github/issues/torresflo/Simple-Math-Library.svg)

<br />
<p align="center">
  <h3 align="center">Simple Math Library</h3>

  <p align="center">
    A simple and light math library that contains usefull classes like Vector, Scalar or Matrix.
    <br />
    <a href="https://github.com/othneildrew/Simple-Math-Library/issues">Report a bug or request a feature</a>
  </p>
</p>

## Table of Contents

* [Getting Started](#getting-started)
  * [Installation](#installation)
* [Usage](#usage)
  * [Example](#example)
* [Contributing](#contributing)
* [License](#license)

## Getting Started

### Installation

Clone the repo (`git clone https:://github.com/torresflo/Simple-Math-Library.git`) and add the Math folder to your project, you are ready to go!

## Usage

### Example

```cpp
#include <Math\Scalar.hpp>
#include <Math\Vector.hpp>
#include <Math\Matrix.hpp>

Math::Vector<float, 4> vector4f; //A vector with 4 floats (x, y, z, w)
Math::Scalar<int, 2> scalar2i; //A scalar with 2 integers
Math::Matrix3x3<double> matrix; //A matrix 3x3 of doubles
```

## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<!-- LICENSE -->
## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.
