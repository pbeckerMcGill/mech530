TODO A5:

APPLYING 3 FAILURE CRITERIA:
(1) Maximum stress
(2) Quadratic Polynomial
(3) HashinFailure Criteria

STEP 1: Load is P1 = 1000N

STEP 2: Solve for on-axis stresses at the top and bottom of each layer

STEP 3: (1) Maximum stress

- Use the 5 Maximum Stress criteria and find the lowest "R" value that will satisfy one of the criteria
- This lowest R-value will be the ratio for which we multiply P1 in order to get first ply failure

STEP 4: (2) Quadratic Polynomial

Calculate the 6 coefficients of the Tsai-Wu Quadratic failure criterion:
- Take on-axis stresses for top and bottom of each layer
- Calculate R * P1
- Sub into R quadratic, solve for R. Take positive root.
- First ply failure is lowest R.

STEP 5: (3) HashinFailure Criteria

- Find R's using this criteria
- Solve for R, take positive values only.
- Failure is lowest R.

STEP 6: PRESENT RESULTS

(1) Maximum stress
- All R values
- Which layer fails first & R value
- Failure mode

(2) Quadratic Polynomial
- All R values
- Which layer fails first & R value
- The load vectors which caused failure

(3) HashinFailure Criteria
- All R values
- Which layer fails first & R value
