# Computing Cox Rings of Embedded Varieties

This repository provides a collection of functions designed to compute Cox rings of subvarieties of Mori dream spaces. These functions implement algorithms to iteratively construct the Cox ring and its grading, leveraging the ambient variety's structure.

---

## Features

### Key Functions

1. **NewGens**  
   - **Description**: Computes the first non-empty B_n of Algorithm 1, which contributes new generators for the Cox ring.  
   - **Input**:  
     - `I`: Ideal of a subvariety.  
     - `lis`: The list m_1,...,m_s.  
     - `h`: Grading map for the generators of the Cox ring of the ambient variety.  
   - **Output**:  
     - B_n \): First non-empty set contributing new variables to the Cox ring.  
     - `lis`: Updated list m_1,...,m_s in the new ring, incorporating new variables for elements of B_n.  
     - `n`: The power adjustment factor for saturation in the algorithm.  
     - `h`: Updated grading map for the generators of the new ring.  

2. **NewRing**  
   - **Description**: Constructs the ideal for a new ring based on input generators and grading.  
   - **Input**:  
     - `I`: Ideal of the subvariety.  
     - `ll`: New generators obtained from the `NewGens` function.  
     - `lis`: List of previous generators m_1,...,m_s.  
     - `n`: Power adjustment factor used in `NewGens`.  
     - `h`: Grading map for the generators of the previous ring.  
   - **Output**:  
     - The ideal of the new ambient ring.  
     - Updated list of generators m_1,...,m_s.  
     - Updated grading map for the new ring.  


3. **IsCox**  
   - **Description**: Checks if the current ideal represents the Cox ring.  
   - **Input**:  
     - `I`: Ideal of the subvariety.  
     - `lis`: The list m_1,...,m_s.
   - **Output**:  
     - Boolean indicating if the ideal is a Cox ring.

4. **ComputeCox**  
   - **Description**: Computes the Cox ring of a variety iteratively.  
   - **Input**:  
     - `I`: Ideal of the subvariety.  
     - `lis`: The list m_1,...,m_s.  
     - `h`: Grading map.  
   - **Output**:  
     - Minimal generators of the Cox ring and their gradings.

5. **IdeCodim2**  
   - **Description**: Computes the codimension 2 irrelevant locus of a toric variety.  
   - **Input**:  
     - `Z`: Toric variety.  
   - **Output**:  
     - Generators of the codimension 2 locus.

6. **DegGens**  
   - **Description**: Computes the degrees of the generators of a given toric variety.  
   - **Input**:  
     - `Z`: Toric variety.  
     - `gens`: List of generators.  
   - **Output**:  
     - Degrees of the generators.

---

## Installation

To use this program:
1. Ensure you have a compatible computational algebra system (e.g., [Magma](https://magma.maths.usyd.edu.au/)).
2. Copy the functions into your workspace or project.
3. Load the functions as required.

---

## Usage

1. **Load Functions**: Include the code in your workspace or project.
2. **Initialize Inputs**: Define the ideal \( I \), the list of generators \( lis \), and the grading map \( h \).
3. **Run Functions**: Call the appropriate function to compute the Cox ring iteratively.

Example usage:
```magma
load "library.m";

Q := GF(32003);  // Define a finite field with 32003 elements.
P := ProjectiveSpace(Q,4);  // Create the 4-dimensional projective space over the finite field Q.
N := Parent(Rays(Fan(P))[1]);  // Extract the lattice associated with the rays of the fan of P.
Z := Blowup(Blowup(P, N![1,1,1,1]), N![-1,0,0,0]);  // Perform two successive blowups on P: first along the ray [1,1,1,1], then along the ray [-1,0,0,0].

QQ := Matrix(Gradings(Z));  // Compute the grading matrix of the Cox ring of the resulting variety Z.
Cl := DivisorClassLattice(Z);  // Obtain the divisor class lattice of Z, representing divisor equivalence classes.
D := Representative(Z, Cl![2,2,3]);  // Define a divisor on Z by choosing a representative from the specified divisor class.

f := &+[Random(Q)*s : s in RiemannRochBasis(D)];  // Construct a random polynomial using a basis of the Riemann-Roch space of the divisor D.
I := Ideal(f);  // Create the ideal defined by the polynomial f.

lis := IdeCodim2(Z);  // Compute the generators of the codimension 2 irrelevant locus of the variety Z.
h := WeilToClassGroupsMap(Z);  // Define the grading map from Weil divisors to the divisor class group of Z.

ComputeCox(I, lis, h);  // Compute the Cox ring of the variety Z using the ideal I, the irrelevant locus lis, and the grading map h.
