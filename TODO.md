# To-Do

## Get project running

- [x] Decide on input data schema
  - May be a good idea to make GSD input in python and then load in C++.  
- [x] Write data parsing scheme
  - [x] Kinematic input
  - [x] Parameter input
- [x] Add logging of startup behavior
- [x] Write data output (preferably in GSD file)
- [x] Figure out why input assertions seem to pass, but data not loaded into class??
- [x] Add logging of startup behavior
- [ ] Write leading order hydrodynamics calculation
  - [ ] Check sign of added mass matrix and gradient
- [ ] Write RKF45 integration function
- [ ] Convert raw pointers into smart pointers

## Additional features

- [ ] Write WCA potential
- [ ] Calculate constraint forces
