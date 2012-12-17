Bicycle Thesis:
===============

**This is for contact forces of constraint bicycle model in steady turning.**

Modules
-------

model.py
---------
Build a Class called bicycle_model which gives benchmark bicycle model with
auxiliary speeds.

steadyturning.py
----------------
Specifical for steady turning maneuver. Here, dynamic equations from forcing
matrix and de_by_inde from nonholonomic equations are key components when 
given a steady-turning configuration of lean and steer angles.

bicycle.py
----------
Here, functions built are only serving for model as well as bicycle parameters
according to Jason' DynamicistToolKit, 
https://github.com/moorepants/DynamicistToolKit.git

steadyturning_conforces.py
---------------------------
Main file in the constraint directory, for the calculation of contact forces
between bicycle wheels and ground.
