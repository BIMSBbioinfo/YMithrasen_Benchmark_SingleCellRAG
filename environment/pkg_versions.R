pk <- readLines("pkgs.txt")
for (p in pk) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(p, as.character(packageVersion(p)), sub(";.*", "", packageDescription(p)$Built), "\n")
  } else cat(p, "NOT INSTALLED\n")
}
