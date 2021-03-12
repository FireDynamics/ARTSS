## Contributing Guidelines

The coding guidelines of ARTSS are mainly based on the [google coding guidelines](https://google.github.io/styleguide/cppguide.html). Accordingly, our rules are not one-to-one with the google standard, therefore we ask you to have a look at the following rules if you want to contribute to ARTSS. It is not necessary to know the google coding guidelines, it is sufficient to read this file.

### Naming Convention
- Variable Naming: We use underscores. The main reason for this is that we have variables which have different meanings depending on upper and lower case (see class Domain.cpp).

  Example:
  
  `int exemplary_variable_name;`

- Naming Branches: Format: `Issue_XY_short_description`

  Example:
  
  `branch: issue_15_coding_guidelines`
  
- Pull Request:
  - Naming Format: `Issue XY short description`

    Example:
  
    `Issue 15 Coding Guidelines`

  - Naming final commit of PR (to master branch): Format: `Issue XY short description (#PR)`. This will be done automatically by GitHub if the PR is named accordingly.

    Example:
  
    `final commit: Issue 15 Coding Guidelines (#113) [description of changes...]`
   
  - Linking to issue. Either link PR to the issue manually or use the phrase "should close #XY" int the description for automatic linking. 

### Parenthesis
- Parenthesis:
  - the opening parenthesis is in the same line as the opening statement and NOT on a separate line
  - do not drop parenthesis even if the the body of a statement consist of a single line

  Example:
  ```
  if (check) {
     body
  }
  ```

### Pointers and References
- Pointers and References: We are using the style where you emphasise the type of the pointed-to data. `someType *somePtr` rather than `someType* somePtr`. The reasoning is mainly that the pointer belongs to the data type because if you define multiple pointers in one line, you have to write: `someType *p1, *p2` while `someType* p1, p2` creates only one pointer.

  Example:
  
  `int *pointer1, *pointer2;`

### Headers
- Self contained: Headers should be self contained, which means it should contain all necessary imports without relying on import of other (imported) headers.

  Example: `Solution.h` uses `std::vector` accordingly it has to import `std::vector`
  
- Include guards: Our include guards have the format ARTSS_DIRNAME_FILENAME_H_

  Example: `Soluction.h` in directory `analysis`
  
  ```
  #ifndef ARTSS_ANALYSIS_SOLUTION_H_
  #define ARTSS_ANALYSIS_SOLUTION_H_
  
  [...]
  
  #endif /* ARTSS_ANALYSIS_SOLUTION_H_ */
  ```
### Line Break/Hard Wrap
 - 100 characters per line is the guideline here but not a rule which has to be followed everywhere. We are adding line breaks if it makes the code more readable, which in most cases means were a breaking a line out of semantic reasons. In which case the line can be much shorter than 100 characters. There are also some special cases where it is not possible to break a line.

As we are only a small group of developers, many parts of the code are not yet compliant with our code guidelines and will be improved accordingly as we make changes. Accordingly, more rules may be added to our coding guidelines if we determine that it is necessary.

Please feel free to contact us if you have any questions.
