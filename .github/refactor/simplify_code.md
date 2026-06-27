# Python Code Simplification Rules


## Files need to be simplied

<!-- - C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\generate_batch_data.py, its task is defined in -->
<!-- C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\.github\tasks\generate_batch_data.md -->
<!-- - C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\centralized_vertical_pca.py, its task is defined in C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\.github\tasks\compute_ver_pca.md -->
<!-- - C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\linear_regression_combinations.py, its task is definded in C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\.github\tasks\compute_lin_reg_products.md -->
<!-- - C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\reconstruction_attack.py and C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\reconstruction_attack.ipynb its task is definded in C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\.github\tasks\reconstruct_ori_data.md -->

- C:\Users\Pham Bao Ngoc\Desktop\7.Semester\Thesis\code\federated_GWAS\plot\plot_solving_times.py, its task is definded in

## Primary Goal

Make the code:

- simpler
- shorter
- easier to read
- easier to maintain

without changing functionality.

## Rules

For additional guidance, follow the rules defined in Section 2 of `copilot-instructions.md`.

### Preserve Behavior

Do not change:

- program output
- business logic
- algorithms
- APIs
- data formats
- tests

unless explicitly requested.

### Prefer Simple Solutions

Use the simplest implementation that correctly solves the problem.

Prefer:

- standard library
- direct logic
- early returns
- small functions
- clear variable names

Avoid:

- unnecessary abstractions
- excessive helper functions
- over-engineering
- premature optimization
- deeply nested code

### Remove Redundancy

Remove:

- unnecessary printing and report
- dead code
- unused imports
- unused variables
- duplicated logic
- redundant comments
- unnecessary wrappers
- unnecessary temporary variables

### Reduce Nesting

Prefer:

```python
if not condition:
    return
