# expts8001.js
Poking about at Matrix indexing mechanisms in Julia

# To do
- performance tests
- different element types
- lower bidiagonal
- sparse
- transpose / adjoint
- update in place vs create new array

# Questions
- indexes method forces a certain iteration order for mult. Does it matter?
- should eg UpperTriangulars be made by inserting into a plain Array and then wrapping at the end?
