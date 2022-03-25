# Generic Core Edge Coupling Proceedure

In the current design we think of the coupling problem as coupling two groups of "solvers". Each solver is either
an internal data computation or as a callback to an external code. For example, three solvers in the GENE/XGC electrostatic
case are 1) GENE compute charge density, 2) XGC compute charge density, 3) XGC field solve for electric field.

Solvers are grouped based on the data transfer needs and a group is typically formed for each high level step in the
coupling. In our concrete example we have the following two groups: 1) compute the charge density and 2) compute the electric field.

The reconciliation process transfers field data from one solver group to another. The field data is received from the first
solver group and combined into a single unified field on the internal intermediate representation, ostensibly Omega_h,
the unified field is converted to, and sent in the native format for the second solver group.


## Basic pseudocode

    // example SolverGroupA is push operation and SolverGroupB is FieldSolve
    Function reconcile(SolverGroupA, SolverGroupB)
        Receive native data for all solvers in SolverGroupA data (core/edge or field-solve)
        FieldTransfer native data for all solvers in SolverGroupA data to Internal(omega_h) mesh
        Construct unified field data on Internal (e.g. take care of buffer/overlap regions)
        FieldTransfer Internal field to SolverGroupB
        Send native data to SolverGroupB
    
    Function solve(SolverGroup)
        asynchronously launch solver for each solver in SolverGroup
        wait for results

    Function run_coupling(Coupling)
        solve(Core/Edge push)
        reconcile(core/edge, field solver)
        solve(field solver)
        // this loop is indicated for case where field solver is operating over multiple partitions and needs to reconcile data between them
        while (solver not converged) 
            reconcile(field solver,field solver)
            solve(field solver)
        reconcile(field solver, core/edge)


    class Coupler {
        Coupler:
            Initialization (receive the mesh data and other initialization data)
        Run:
            while (not done coupling)
                // asyncrounously run each of the coupling steps
                run_coupling(charge density+electric field coupling) // poisson
                run_coupling(current density+parallel part of vector potential coupling) // ampere
                ...
                // wait for all of the coupling steps to complete
    }

### Functions
1. solve: launch the solver for the Core/Edge coupling this corresponds to either the push operation, or the field solve.
2. reconcile: Gets data from one solver's/coordinate system's and converts it to another's coordinate systems by way of intermediate representation
3. run_coupling: Takes solve and reconcile policy and performs the main control flow

### Classes/Concepts
1. CoordinateSystem
2. Field
3. Solver
4. FieldTransfer
5. DataTransfer : how data transfer will happen for a given code
6. SolverGroup : group of solvers
7. Coupling : defines a complete coupling procedure for a given field/group of solvers e.g. this could be the charge density/poisson solve
8. Coupler : Ties together couplings, data initialization, etc to provide the full

### open questions
1. For coupling with multiple fields should solve,reconcile, etc. deal with multiple fields internally? Or, run_coupling
    is run independently for each coupled set of fields? I think the cleaner approach is the latter where we run the
    coupling for each set of fields. If it becomes necessary later we can push the fields down the stack. Only problem
    I foresee is if there is some shared metadata or state. This may make parallel easier as we can launch the run_coupling
    for each coupled field on a new thread or process.
2. Design for async function calls. Clearly there is a lot of the design that must be done async and execute in parallel
    for any hope of efficiency. Question is how do we want to approach that? We can just dump stuff onto raw threads, 
    look at sender/receiver impl. in HPX/reference impl, coroutines, futures, etc. 
3. Who owns the state of the state of when the coupling is complete?

## Notes and things to consider
1. Couplings may need to be performed at different rates.
2. Since couplings are performed asynchronously, then they must be orthogonal in the
    shared data they write to avoid race conditions.
3. Each solver in the solver group needs to also have a field transfer method to go to
   the internal representation. Needs to be part of the solver because each solver may
   have a different way of getting to the internal state.