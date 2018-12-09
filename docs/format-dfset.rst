How to make a DFSET file
------------------------

.. _label_format_DFILE:

Format of ``DFSET``
~~~~~~~~~~~~~~~~~~~

The displacement-force data sets obtained by first-principles (or classical force-field) calculations
have to be saved to a file, say *DFSET*. Then, the force constants are estimated by setting ``DFSET = `` *DFSET* and with ``MODE = optimize``.

The *DFSET* file must contain the atomic displacements and corresponding forces in Cartesian coordinate for at least ``NDATA`` structures (displacement patterns)
in the following format: 

.. math::
    :nowrap:
     
    \begin{verbatim}
      # Structure number 1 (this is just a comment line)
    \end{verbatim}
    \begin{eqnarray*}
     u_{x}(1) & u_{y}(1) & u_{z}(1) & f_{x}(1) & f_{y}(1) & f_{z}(1) \\
     u_{x}(2) & u_{y}(2) & u_{z}(2) & f_{x}(2) & f_{y}(2) & f_{z}(2) \\
              & \vdots   &          &          & \vdots   &          \\
     u_{x}(\mathrm{NAT}) & u_{y}(\mathrm{NAT}) & u_{z}(\mathrm{NAT}) & f_{x}(\mathrm{NAT}) & f_{y}(\mathrm{NAT}) & f_{z}(\mathrm{NAT})
    \end{eqnarray*}
    \begin{verbatim}
      # Structure number 2 
    \end{verbatim}
     \begin{eqnarray*}
     u_{x}(1) & u_{y}(1) & u_{z}(1) & f_{x}(1) & f_{y}(1) & f_{z}(1) \\
       & \vdots   &          &          & \vdots   &          \\
     \end{eqnarray*}


When there are ``NAT`` atoms in the supercell and ``NDATA`` data sets, 
there should be  ``NAT`` :math:`\times` ``NDATA`` lines in the ``DFSET``.
The unit of displacements and forces must be **Bohr** and **Ryd/Bohr**, respectively.


Generation of ``DFSET`` by extract.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

