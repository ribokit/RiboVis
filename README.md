# pymol_daslab

![1q9a example image](https://raw.github.com/DasLab/pymol_daslab/master/1q9a.png)

**pymol_daslab** is a set of useful and short *Python* functions for making pictures of RNA and proteins in pymol in our 'lab style'. For a quick preview:

| Function | Description |
| --- | --- |
| `rr()` | For RNA: with 2&prime; OH as spheres, bases as filled rings, and backbone as cartoon ribbons, rainbow colored from 5&prime; to 3&prime;. No hydrogens, white background. |
| `rd()` | For proteins: side chains are all-atom and colored _CPK_, backbone is rainbow cartoon from N to C terminus. |
| `sa()` | Superimposes all models to the first one. _Thanks to Kyle Beauchamp for this one_ |

And more! ...

## Installation

To install **pymol_daslab**, simply:

- From GitHub, download the zip or tar file of the repository and unpack; or 

```bash
git clone https://github.com/DasLab/pymol_daslab.git
```

- In _PyMol_, type:

```python
run pymol_daslab.py
```

Tested quickly with:

```python
fetch 1q9a
rr()
```

- _(Optional)_ Create or edit a `.pymolrc` file in your home directory, add these lines:

```python
import sys
sys.path.append('/path/to/pymol_daslab')
run /path/to/pymol_daslab/pymol_daslab.py
```

> Replce with your `/path/to/pymol_daslab`.

This will automatically load `pymol_daslab` upon start every time.

## Documentation

Documentation is available at https://ribokit.github.io/pymol_daslab/.

## License

Copyright &copy; of **pymol_daslab** _Source Code_ is described in [LICENSE.md](https://github.com/DasLab/pymol_daslab/blob/master/LICENSE.md).

<br/>
Developed by **Das lab**, _Leland Stanford Junior University_.
<br/>
README by [**t47**](http://t47.io/), *April 2016*.

