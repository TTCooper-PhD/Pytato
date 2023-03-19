# Available Enzymes in Pyteomics

The [Pyteomics library](https://github.com/levitsky/pyteomics) provides support for several proteolytic enzymes. Here is a list of enzymes available in the library, along with their case-sensitive names and cleavage rules:

1. **trypsin**: `(?<=[KR])(?!P)`
   - Cleaves after lysine (K) and arginine (R) residues, unless followed by a proline (P).

2. **lys-c**: `(?<=K)(?!P)`
   - Cleaves after lysine (K) residues, unless followed by a proline (P).

3. **lys-n**: `(?=K)`
   - Cleaves before lysine (K) residues.

4. **arg-c**: `(?<=R)(?!P)`
   - Cleaves after arginine (R) residues, unless followed by a proline (P).

5. **asp-n**: `(?=D)`
   - Cleaves before aspartic acid (D) residues.

6. **glu-c**: `(?<=[DE])(?!P)`
   - Cleaves after glutamic acid (E) and aspartic acid (D) residues, unless followed by a proline (P).

7. **chymotrypsin**: `(?<=[FWY])(?!P)`
   - Cleaves after phenylalanine (F), tryptophan (W), and tyrosine (Y) residues, unless followed by a proline (P).

8. **thermolysin**: `(?<=[LIFVAM])`
   - Cleaves after leucine (L), isoleucine (I), phenylalanine (F), valine (V), alanine (A), and methionine (M) residues.

Remember that the enzyme names are case-sensitive when used with the Pyteomics parser.
