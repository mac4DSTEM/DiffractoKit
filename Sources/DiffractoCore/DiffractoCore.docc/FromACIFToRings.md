# From a CIF to rings

Read a structure from CIF text, see the symmetry the kit finds in it, and list the reflections and rings below a
reciprocal radius.

## Read the structure

``CIFImport/read(_:fileBaseName:)`` takes the text of a CIF and returns a ``CIFDocument``; its `structure` is a ``CrystalStructure``.
A file the kit cannot trust throws a ``CIFImportError`` whose description names the tag or value at fault.

```swift
import DiffractoCore

let gold: CrystalStructure = try CIFImport.read(cifText).structure
let laue: String = gold.laueClass?.symbol ?? "?"
print(gold.formula, laue, gold.operations.count)      // Au4 m-3m 192
let dFirst: Double = gold.cell.d(MillerIndex(1, 1, 1))   // 2.3545 A
```

The symmetry is found from the atoms, not read from the file: the file's own operations are kept as `fileOperations`, and
`fileOperationsAreSymmetries` says whether they were right. A structure can also be built in code, from
``UnitCell/cubic(a:)`` and a list of ``AtomSite``.

## List the reflections

``CrystalStructure/reflections(kMax:)`` returns a ``ReflectionList``: every hkl with 0 < |g| <= k_max, in a fixed order,
each flagged by ``Reflection/isAllowed`` when the structure's operations do not forbid it. Lengths are in angstroms and
k = 1/d is in inverse angstroms, so `kMax: 1.2` reaches d = 0.83 A.

```swift
let list: ReflectionList = try gold.reflections(kMax: 1.2)
for family in list.families(allowedOnly: true).prefix(4) {
    print(family.name, family.d, family.multiplicity)   // {1 1 1} 2.3545 8, {2 0 0} 2.0391 6, ...
}
```

A ``ReflectionFamily`` is one ring: the reflections whose |g| agree to ``ReflectionList/familyTolerance``, named by the
orbit representatives (``MillerIndex`` ordering), with ``ReflectionFamily/d`` and ``ReflectionFamily/multiplicity``.
A request that would enumerate too much (a k_max above ``ReflectionList/maximumKMax``, or a cell so large that the list
passes ``ReflectionList/maximumReflections``) throws a ``ReflectionListError`` instead of returning a partial list.

## Next

Allowed is not the same as visible: the intensity of each ring comes from the structure factors in DiffractoScattering,
which cuts the list again by |F|.
