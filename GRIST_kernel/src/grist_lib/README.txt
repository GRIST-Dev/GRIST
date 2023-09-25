How to separate "model" from "infrastructure" (some thoughts)
1. Features of infrastructure include:
(i) definition of data strucures and software functions (e.g., domain decomposision, i/o and communication).
(ii) shared math and well-defined rules (e.g., time manager, error handling).
(iii) common namelist for different models.
A basic rule to think is that whether the infrastructure can support more than one model.

2. Features of "model" include:
(i) decribe a certain type of target (e.g., 1d profile atmosphere/scm; 2d atmosphere/ocean/shallow water; 3d atmosphere/ocean).
(ii) define a structral relationship of various components to yield a model (e.g., dycore-tracer-physics).
(iii) within each model, the computational logic, algorithms are left to modelers' ideas to descide.

Some thinking can be motivated from an interesting discussion related to modularity:
Held, I., and D. Randall, (2011), Point/Counterpoint. IEEE Software, 28(6), 62-65.
