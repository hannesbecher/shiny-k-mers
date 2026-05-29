# AGENTS.md file for Tetmer

## Documenting instuctions
- use Roxygen2 markdown comments to document functions and classes in the codebase
- Do not update the man directory, as it is generated automatically from the Roxygen2 comments in the codebase
- Do not update the NAMESPACE file, as it is also generated automatically from the Roxygen2 comments in the codebase
- do not add method registrations by hand, handle these with Roxygen2 comments in the codebase
- run roxygen2::roxygenise(clean=TRUE) to generate docs