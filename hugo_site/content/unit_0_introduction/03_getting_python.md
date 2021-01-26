---
title: "Getting Python"
date: 2020-11-24T16:52:22Z
draft: false
pre: "3. "
---

{{% notice warning %}}
This course and all of the sample solutions use Python 3.
{{% /notice %}}

## Windows

You can download and install Python 3 using the releases provided by the CPython team 
[here](https://www.python.org/downloads/windows/). Installation using the windows 
installer packages should be straightforward, but if required you can find further 
documentation on the installation process 
[here](https://docs.python.org/3/using/windows.html).

## Mac OS X

You can download and install Python 3 on a Mac using the Python [download 
page](https://www.python.org/downloads/mac-osx/), and there are further instructions in 
the [official documentation](https://docs.python.org/3/using/mac.html).

As an alternative, Homebrew is a package manager for Mac OSX that is very useful for 
installing scientific software and dependencies. You can see the official Homebrew 
[page](https://brew.sh/) for installation instructions. Once Homebrew is installed you 
can install Python using

```bash
brew install python
```

## Debian and Ubuntu 20.04 Linux

Debian Linux, and Ubuntu 20.04 both ship with Python 3 pre-installed. 


## Virtual environments

As with everything you do in Python, it is recommended that you work in a virtual 
environment during this course. You can create a virtual environment in Python 3 using 
the `venv` module

```bash
python3 -m venv env
```

This creates a directory `env` containing the virtual environment. You can then activate 
this environment using

```bash
source env/bin/activate
```

Once your environment is activated, you can install the dependencies for this course 
using `pip`

```bash
pip install numpy scipy matplotlib
```






