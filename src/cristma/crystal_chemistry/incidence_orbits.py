"""Oriented local contact incidences derived from undirected pair orbits."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import math

from cristma.crystallography import (
    AsymmetricUnitMapping,
    PeriodicSymmetryRelation,
    SymmetryContext,
    SymmetryPairTable,
    periodic_endpoint_instance,
)
from cristma.crystallography.symmetry_context import _digest
from cristma.structure import CrystalStructure

from .contacts import ResolutionStatus, SecondaryEvidence
from .orbit_contacts import (
    ContactInterpretation,
    EndpointRole,
    OrientationMode,
    ResolvedContactOrbit,
)


@dataclass(frozen=True, slots=True)
class ContactIncidenceOrbit:
    incidence_orbit_id: str
    resolved_contact_orbit_id: str
    interpretation_id: str
    center_independent_site_id: str
    ligand_independent_site_id: str
    oriented_periodic_relation: PeriodicSymmetryRelation
    equivalent_oriented_relations: tuple[PeriodicSymmetryRelation, ...]
    incidence_multiplicity_per_center: int
    effective_neighbor_occupancy: float
    status: ResolutionStatus
    evidence: tuple[SecondaryEvidence, ...]

    def __post_init__(self) -> None:
        identities = (
            self.incidence_orbit_id,
            self.resolved_contact_orbit_id,
            self.interpretation_id,
            self.center_independent_site_id,
            self.ligand_independent_site_id,
        )
        if any(not value for value in identities):
            raise ValueError("contact incidence identities must not be empty")
        if self.incidence_multiplicity_per_center <= 0:
            raise ValueError("incidence multiplicity must be positive")
        if (
            tuple(sorted(set(self.equivalent_oriented_relations)))
            != self.equivalent_oriented_relations
            or not self.equivalent_oriented_relations
        ):
            raise ValueError("oriented incidence relations must be non-empty, unique, and sorted")
        if self.oriented_periodic_relation != self.equivalent_oriented_relations[0]:
            raise ValueError("incidence representative must be its first oriented relation")
        if self.incidence_multiplicity_per_center != len(self.equivalent_oriented_relations):
            raise ValueError("incidence multiplicity must equal its exact relation count")
        if not math.isfinite(self.effective_neighbor_occupancy) or not (
            0.0 <= self.effective_neighbor_occupancy <= 1.0 + 1e-12
        ):
            raise ValueError("effective neighbor occupancy must lie between zero and one")


@dataclass(frozen=True, slots=True)
class _OrientedSeed:
    center_site_id: str
    ligand_site_id: str
    relation: PeriodicSymmetryRelation
    ligand_component_side: int


def _seeds(
    geometry,
    interpretation: ContactInterpretation,
    context,
) -> tuple[_OrientedSeed, ...]:
    first = geometry.first_independent_site_id
    second = geometry.second_independent_site_id
    relation = geometry.canonical_relation
    if interpretation.orientation_mode is OrientationMode.UNDIRECTED:
        return (
            _OrientedSeed(first, second, relation, 1),
            _OrientedSeed(second, first, relation.inverse(context), 0),
        )
    if interpretation.endpoint_roles == (EndpointRole.CENTER, EndpointRole.LIGAND):
        return (_OrientedSeed(first, second, relation, 1),)
    if interpretation.endpoint_roles == (EndpointRole.LIGAND, EndpointRole.CENTER):
        return (_OrientedSeed(second, first, relation.inverse(context), 0),)
    raise ValueError("oriented contact interpretation has invalid endpoint roles")


def _effective_ligand_occupancy(
    interpretation: ContactInterpretation,
    seeds: tuple[_OrientedSeed, ...],
    structure: CrystalStructure,
) -> float:
    participating_species = set()
    for seed in seeds:
        for record in interpretation.component_pair_interpretations:
            participating_species.add(
                record.second_species if seed.ligand_component_side == 1 else record.first_species
            )
    ligand_id = seeds[0].ligand_site_id
    ligand_site = next(site for site in structure.sites if site.id == ligand_id)
    return math.fsum(
        float(component.occupancy.value)
        for component in ligand_site.components
        if component.species in participating_species
    )


class ContactIncidenceBuilder:
    """Build local incidence identities without selecting shell boundaries."""

    def build(
        self,
        pair_table: SymmetryPairTable,
        contact_orbits: tuple[ResolvedContactOrbit, ...],
        structure: CrystalStructure,
        mapping: AsymmetricUnitMapping,
        context: SymmetryContext,
    ) -> tuple[ContactIncidenceOrbit, ...]:
        if not isinstance(pair_table, SymmetryPairTable):
            raise TypeError("pair_table must be SymmetryPairTable")
        if not isinstance(structure, CrystalStructure):
            raise TypeError("structure must be CrystalStructure")
        if not isinstance(mapping, AsymmetricUnitMapping):
            raise TypeError("mapping must be AsymmetricUnitMapping")
        if not isinstance(context, SymmetryContext):
            raise TypeError("context must be SymmetryContext")
        if pair_table.asymmetric_unit_mapping_fingerprint != mapping.fingerprint:
            raise ValueError("pair table and asymmetric-unit mapping disagree")
        if pair_table.symmetry_context_fingerprint != context.fingerprint:
            raise ValueError("pair table and symmetry context disagree")
        geometries = {
            orbit.geometry_orbit_id: orbit for orbit in pair_table.contact_orbits
        }
        output: list[ContactIncidenceOrbit] = []
        for resolved in contact_orbits:
            try:
                geometry = geometries[resolved.geometry_orbit_id]
            except KeyError as exc:
                raise ValueError("resolved contact references an unknown geometry orbit") from exc
            for interpretation in resolved.interpretations:
                grouped: dict[tuple[str, str], list[_OrientedSeed]] = defaultdict(list)
                for seed in _seeds(geometry, interpretation, context):
                    grouped[(seed.center_site_id, seed.ligand_site_id)].append(seed)
                for (center_id, ligand_id), raw_seeds in grouped.items():
                    seed_group = tuple(raw_seeds)
                    neighbor_relations: dict[
                        tuple[str, str, tuple[int, int, int]],
                        PeriodicSymmetryRelation,
                    ] = {}
                    center_stabilizer = mapping.by_site_id[center_id].stabilizer_relations
                    ligand_stabilizer = mapping.by_site_id[ligand_id].stabilizer_relations
                    for seed in seed_group:
                        for center_relation in center_stabilizer:
                            moved = center_relation.compose(seed.relation, context)
                            canonical = min(
                                moved.compose(ligand_relation, context)
                                for ligand_relation in ligand_stabilizer
                            )
                            endpoint = periodic_endpoint_instance(
                                ligand_id,
                                canonical,
                                mapping,
                            )
                            previous = neighbor_relations.get(endpoint)
                            if previous is None or canonical < previous:
                                neighbor_relations[endpoint] = canonical
                    relations = tuple(sorted(neighbor_relations.values()))
                    if not relations:
                        raise ValueError("contact incidence has no local neighbor relation")
                    incidence_id = "contact-incidence-orbit:" + _digest(
                        {
                            "resolved_contact_orbit_id": resolved.resolved_contact_orbit_id,
                            "interpretation_id": interpretation.interpretation_id,
                            "center_site_id": center_id,
                            "ligand_site_id": ligand_id,
                            "relations": tuple(
                                (item.operation_key, item.lattice_translation)
                                for item in relations
                            ),
                        }
                    )
                    status = (
                        ResolutionStatus.INCOMPLETE
                        if resolved.status is ResolutionStatus.INCOMPLETE
                        or interpretation.status is ResolutionStatus.INCOMPLETE
                        else ResolutionStatus.AMBIGUOUS
                        if resolved.status is ResolutionStatus.AMBIGUOUS
                        else ResolutionStatus.RESOLVED
                    )
                    output.append(
                        ContactIncidenceOrbit(
                            incidence_id,
                            resolved.resolved_contact_orbit_id,
                            interpretation.interpretation_id,
                            center_id,
                            ligand_id,
                            relations[0],
                            relations,
                            len(relations),
                            _effective_ligand_occupancy(
                                interpretation,
                                seed_group,
                                structure,
                            ),
                            status,
                            interpretation.evidence,
                        )
                    )
        return tuple(sorted(output, key=lambda item: item.incidence_orbit_id))
__all__ = ["ContactIncidenceBuilder", "ContactIncidenceOrbit"]
