# expts8001.jl
Poking about at Matrix indexing mechanisms in Julia, looking for
solutions to https://github.com/JuliaLang/julia/issues/8001

# To do
- performance tests
- different element types
- DONE lower bidiagonal, tridiagonal
- DONE convert fully to using Banded to represent structure
  But backed out for using structure to determine target matrix type
  because that appears to defeat the inliner
- experiment with copying/reshaping original array(s) to make indexing faster
- sparse
- transpose / adjoint
- reduce cost of indexing second array in dot() - make sure inline and not bounds checked
- update in place vs create new array
- subarrays
- look into generated functions

# Questions
- indexes method forces a certain iteration order for mult. Does it matter?
- should eg UpperTriangulars be made by inserting into a plain Array and then wrapping at the end?
