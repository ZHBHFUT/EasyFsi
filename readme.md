<a name="readme-top"></a>
<!--
*** This readme is written based on Best-README-Template.
-->

<!-- PROJECT SHIELDS -->

<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/ZHBHFUT/EasyFsi">
    <img src="images/logo.png" alt="Logo" width="155" height="144">
  </a>

  <h3 align="center">EasyFsi Library</h3>

  <p align="center">
    An flexable toolkit for solving FSI or multi-physics coupling problem!
    <br />
    <a href="https://github.com/ZHBHFUT/EasyFsi/blob/master/EasyFsi.md"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/ZHBHFUT/EasyFsi/blob/master/demo.cpp">View Demo</a>
    <a href="https://github.com/ZHBHFUT/EasyFsi/issues">Report Bug</a>
    <a href="https://github.com/ZHBHFUT/EasyFsi/issues">Request Feature</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

This library provides a flexible toolkit for integrating different solvers to solve fluid-structure interaction and other multiphysics coupling problems. It can be used to solve aeroelastic, aerothermoelastic, etc.

Why use this library?
+ Suitable for MPI parallel computing, which automatically assembles distributed boundary
+ Implement many interpolation methods such as spline interpolation, ISO-parameteric inversion, geometric mapping, etc.
+ Communication between different solvers based on socket
+ No dependencies except pybind11
+ Easily integrated into C/C++, Fortran, Python, MATLAB solvers.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Built With

* [Visual Studio][visualstudio-url]

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

### Prerequisites

* Visual Studio, make sure C++20 is supported
* pybind11, if you need python binding

### Installation

_Below is an example of how you can instruct your audience on installing and setting up your app. This template doesn't rely on any external dependencies or services._

1. Clone the repo
   ```sh
   git clone https://github.com/ZHBHFUT/EasyFsi.git
   ```
2. Open `EasyFsi.sln` and build all projects.
3. Include `EasyFsi.h` to your C/C++ project, see `demo.cpp` for usage.
4. Use `import EasyFsi` in your python script, see `demo.py` for usage.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

See [`demo.cpp`](https://github.com/ZHBHFUT/EasyFsi/blob/master/demo.cpp) and [`demo.py`](https://github.com/ZHBHFUT/EasyFsi/blob/master/demo.cpp) to find out how to use this library.

More documents please see [`EasyFsi.md`](https://github.com/ZHBHFUT/EasyFsi/blob/master/EasyFsi.md).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ROADMAP -->
## Roadmap

- [ ] Add Changelog
- [ ] Add MATLAB binding
- [ ] Implement the `CoupledRegion` class
- [ ] Add additional interpolation method
- [ ] Add document

See the [open issues](https://github.com/ZHBHFUT/EasyFsi/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

ZHANG Bing - zhangbing@hfut.edu.cn

Project Link: [https://github.com/ZHBHFUT/EasyFsi](https://github.com/ZHBHFUT/EasyFsi)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

TBD

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[license-shield]: https://github.com/ZHBHFUT/EasyFsi/blob/master/images/LICENSE-MIT.svg?style=for-the-badge
[license-url]: https://mit-license.org/
[visualstudio-img]: https://visualstudio.microsoft.com/wp-content/uploads/2021/10/Product-Icon.svg
[visualstudio-url]: https://visualstudio.microsoft.com
[pybind11-img]:https://github.com/pybind/pybind11/raw/master/docs/pybind11-logo.png
[pybind11-url]:https://github.com/pybind/pybind11
