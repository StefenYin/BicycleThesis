==============
Bicycle Thesis
==============

**This is for contact forces of constraint bicycle model in steady turning.**

Modules
=======

model
-----
Build a Class called **bicycle_model** which gives benchmark bicycle model with
auxiliary speeds.

steadyturning
-------------
Specifical for steady turning maneuver. Here, *dynamic_equations* from forcing
matrix and *de_by_inde* from nonholonomic equations are key components when 
given a steady-turning configuration of lean and steer angles.

bicycle
-------
Here, functions built are only serving for model as well as bicycle parameters
according to Jason' **DynamicistToolKit**, 
https://github.com/moorepants/DynamicistToolKit.git

conforces_steadyturning
-----------------------
Main file in the **constraint** directory, for the calculation of contact forces
between bicycle wheels and ground. A Class named **steady_turning** is built.

test_model, test_steadyturning, test_bicycle
--------------------------------------------
Test modules, but did not use assertation yet.


Folders
=======

seperated_file
--------------
Containing original model, and contact forces files.
