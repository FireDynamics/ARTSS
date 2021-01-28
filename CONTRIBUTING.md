## Contributing Guidelines

The coding guidelines of ARTSS are mainly based on the [google coding guidelines](https://google.github.io/styleguide/cppguide.html). Below there are some rules which are not defined in the google standard or we want to especially emphasise them in order to achieve the most consistent code possible.

- Parenthesis:
  - the opening parenthesis is in the same line as the opening statement and NOT on a separate line
  - do not drop parenthesis even if the the body of a statement consist of a single line

Example:
```
if (check) {
   body
}
```

- Naming Convention: We use underscores. The main reason for this is that we have variables which have different meanings depending on upper and lower case (see class Domain.cpp).

- Pointers and References: We are using the style where you emphasise the type of the pointed-to data. `someType *somePtr` rather than `someType* somePtr`. The reasoning is mainly that the pointer belongs to the data type because if you define multiple pointers in one line, you have to write: `someType *p1, *p2` while `someType* p1, p2` creates only one pointer.


As we are only a small group of developers, many parts of the code are not yet compliant with our code guidelines and will be improved accordingly as we make changes. Accordingly, more rules may be added to our coding guidelines if we determine that it is necessary.

Please feel free to contact us if you have any questions.
