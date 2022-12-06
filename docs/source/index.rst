Welcome to Miritools's documentation!
===================================

The **miritools** package contains convenience functions usefull to display or characterize MIRI data on board JWST. Some are more generic and could be usefull in general, but all of them have been created with MIRI in mind in some way or another.

Miritools start at v4.0.0. This is for historical reasons because it is based on the miricap package and I kept the version number for convenience.

Check out the :doc:`sub-package` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   sub-package
   api


.. _installation:

How to install
------------------
The latest version of this software can be retrieved from the repository hosted at:
https://github.com/ccossou/miritools

This package is a simple Python module, and all default commands
works, such as:

..  code-block:: console

    pip install miritools

How to use
-----------------
* First import the package:

..  code-block:: python

    import miritools
    # or
    from miritools import mask, read, plot, flux  #...

* To get help on a given function, either look at the documentation or use the built-in *help* function:

..  code-block:: python

    from miritools import plot
    help(plot.single_image)

