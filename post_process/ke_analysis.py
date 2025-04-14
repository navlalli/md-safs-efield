""" Analyze kinetic energies """

import numpy as np
import math
import os

def overview(dump_file):
    """ Find number of timesteps and domain properties in dump file """
    # Find number of lines in file
    with open(dump_file) as f:
        for count, _ in enumerate(f):
            pass
        nlines = count + 1
    # Find number of particles and number of timesteps in file
    with open(dump_file) as f:
        for _ in range(3):
            f.readline()
        natoms = int(f.readline())
        ntimes = round(nlines / (natoms + 9))
        f.readline()
        lsplit = f.readline().split()
        xlo, xhi = float(lsplit[0]), float(lsplit[1])
        lsplit = f.readline().split()
        ylo, yhi = float(lsplit[0]), float(lsplit[1])
        lsplit = f.readline().split()
        zlo, zhi = float(lsplit[0]), float(lsplit[1])
    assert math.isclose(xlo, ylo)
    assert math.isclose(xlo, zlo)
    assert math.isclose(xhi, yhi)
    assert math.isclose(xhi, zhi)
    box_len = xhi - xlo  # Cubic box length (Angstrom)
    return ntimes, natoms, box_len

def mass_array(dump_file, natoms):
    """ Create a mass array for vectorised calculations """
    mass = {}
    mass["H"] = 1.0080
    mass["C"] = 12.0107
    mass["O"] = 15.9994
    mass_arr = np.zeros(natoms)
    with open(dump_file) as f:
        for _ in range(9):
            f.readline()
        for i in range(natoms):
            lsplit = f.readline().split()
            element = lsplit[1]
            mass_arr[i] = mass[element]
    return mass_arr

def find_molecules(mol_ids):
    """ Analyze mol ids to find molecules """
    assert mol_ids[0] == 1
    mol_ids_unique, mol_atom_counts = np.unique(mol_ids, return_counts=True)
    nmolecules = len(mol_ids_unique)
    assert len(mol_atom_counts) == nmolecules
    assert np.max(mol_ids) == nmolecules
    return mol_ids_unique, mol_atom_counts, nmolecules
    
def rot_analysis_two_atoms(angmom, inertia, xyzsplit, check=0):
    """ Find rotational kinetic energy and angular velocity for molecules consisting
    of two atoms 
    """
    # Select three principal axes using Gram-Schmidt process
    # e1 is the axis with 0 principal moment of inertia
    e1 = xyzsplit[1, :] - xyzsplit[0, :]
    e1 /= np.linalg.norm(e1)
    e2 = np.random.randn(3)
    e2 -= e2.dot(e1) * e1
    e2 /= np.linalg.norm(e2)
    e3 = np.cross(e1, e2)
    eig_vec = np.column_stack((e1, e2, e3))
    # Find principal moments of inertia - I1 = 0.0
    eig = np.zeros(3)
    eig_tmp = np.matmul(inertia, e2) / e2 
    assert np.all(np.isclose(eig_tmp, eig_tmp[0]))
    eig[1] = eig_tmp[0]
    eig_tmp = np.matmul(inertia, e3) / e3 
    assert np.all(np.isclose(eig_tmp, eig_tmp[0]))
    eig[2] = eig_tmp[0]

    # Find rotational KE and angular velocity of molecule: w1 = 0.0
    rot_ke = 0.0
    ang_vel = np.zeros(3)
    pr_ang_vel = np.zeros(3)
    for i in range(1, 3):
        rot_ke += 0.5 * np.sum(angmom * eig_vec[:, i])**2 / eig[i]
        pr_ang_vel[i] += np.sum(angmom * eig_vec[:, i]) / eig[i]
    for i in range(3):
        ang_vel[0] += pr_ang_vel[i] * eig_vec[0, i]
        ang_vel[1] += pr_ang_vel[i] * eig_vec[1, i]
        ang_vel[2] += pr_ang_vel[i] * eig_vec[2, i]

    if check:
        tol = 1e-6
        eig, eig_vec = np.linalg.eigh(inertia)  # Principal moments of inertia and directions
        assert np.min(eig) > -tol, f"{np.min(eig) = }"
        max_eig = np.max(eig)
        thres = tol * max_eig  # Principal angular velocity is 0 when I < thres
        rot_ke_test = 0.0
        pr_ang_vel_test = np.zeros(3)
        ang_vel_test = np.zeros(3)
        for i in range(3):
            if eig[i] > thres:
                rot_ke_test += 0.5 * np.sum(angmom * eig_vec[:, i])**2 / eig[i]
                pr_ang_vel_test[i] += np.sum(angmom * eig_vec[:, i]) / eig[i]
        assert math.isclose(rot_ke_test, rot_ke, abs_tol=1e-7), f"{rot_ke_test = }, {rot_ke = }"
        for i in range(3):
            ang_vel_test[0] += pr_ang_vel_test[i] * eig_vec[0, i]
            ang_vel_test[1] += pr_ang_vel_test[i] * eig_vec[1, i]
            ang_vel_test[2] += pr_ang_vel_test[i] * eig_vec[2, i]
        assert np.all(np.isclose(ang_vel, ang_vel_test))
        rot_ke_test2 = 0.5 * np.sum(angmom * ang_vel_test)
        assert math.isclose(rot_ke_test2, rot_ke, abs_tol=1e-6)

    return rot_ke, ang_vel

def rot_analysis(angmom, inertia, check=0):
    """ Find rotational kinetic energy from angular momentum and moment of inertia 
    tensor - special treatment for linear molecules.
    """
    tol = 1e-6
    if np.linalg.det(inertia) > tol:
        ang_vel = np.matmul(np.linalg.inv(inertia), np.transpose(angmom))
        rot_ke = 0.5 * np.sum(ang_vel * angmom)
        if check:
            eig, eig_vec = np.linalg.eigh(inertia)  # Principal moments and axes
            rot_ke_test = 0.0
            pr_ang_vel = np.zeros(3)  # Principal angular velocities
            ang_vel_test = np.zeros(3)  # Principal angular velocities
            for i in range(3):
                rot_ke_test += 0.5 * np.sum(angmom * eig_vec[:, i])**2 / eig[i]
                pr_ang_vel[i] += np.sum(angmom * eig_vec[:, i]) / eig[i]
            for i in range(3):
                ang_vel_test[0] += pr_ang_vel[i] * eig_vec[0, i]
                ang_vel_test[1] += pr_ang_vel[i] * eig_vec[1, i]
                ang_vel_test[2] += pr_ang_vel[i] * eig_vec[2, i]
            assert math.isclose(rot_ke, rot_ke_test, abs_tol=1e-6)
            assert math.isclose(ang_vel[0], ang_vel_test[0], abs_tol=1e-6)
            assert math.isclose(ang_vel[1], ang_vel_test[1], abs_tol=1e-6)
            assert math.isclose(ang_vel[2], ang_vel_test[2], abs_tol=1e-6)
    else:  # Perhaps for a longer linear molecule, such as CO2
        eig, eig_vec = np.linalg.eigh(inertia)  # Principal moments of inertia and directions
        assert np.min(eig) > -tol, f"{np.min(eig) = }"
        max_eig = np.max(eig)
        thres = tol * max_eig  # Principal angular velocity is 0 when I < thres
        rot_ke = 0.0
        pr_ang_vel = np.zeros(3)
        ang_vel = np.zeros(3)
        for i in range(3):
            if eig[i] > thres:
                rot_ke += 0.5 * np.sum(angmom * eig_vec[:, i])**2 / eig[i]
                pr_ang_vel[i] += np.sum(angmom * eig_vec[:, i]) / eig[i]
        for i in range(3):
            ang_vel[0] += pr_ang_vel[i] * eig_vec[0, i]
            ang_vel[1] += pr_ang_vel[i] * eig_vec[1, i]
            ang_vel[2] += pr_ang_vel[i] * eig_vec[2, i]
        if check:
            rot_ke_test = 0.5 * np.sum(angmom * ang_vel)
            assert math.isclose(rot_ke_test, rot_ke, abs_tol=1e-6)
    return rot_ke, ang_vel

def mol_analysis(mol_ids, element, mass_arr, xyz_uw, vxyz):
    """ Analyze molecules """
    # Find molecules
    mol_ids_unique, mol_atom_counts, nmolecules = find_molecules(mol_ids)
    # Analyze velocities
    com = np.zeros((nmolecules, 3))  # Position of centre of mass of each molecule (Ang)
    vcom = np.zeros((nmolecules, 3))  # Velocity of centre of mass of each molecule (Ang/fs)
    ang_vel = np.zeros((nmolecules, 3))  # Angular velocity of each molecule (/fs)
    type_counts = np.zeros((nmolecules, 3), dtype=int)  # Counts of H, C, O in each molecule
    mass_mol = np.zeros(nmolecules)  # Mass of each molecule (g/mol)
    tot_ke_mol = np.zeros(nmolecules)  # Total KE of each molecule
    tot_ke_mol_cpts = np.zeros((nmolecules, 3))  # Total KE of each molecule - xyz cpts
    rot_ke_mol = np.zeros(nmolecules)  # Rotational KE of each molecule
    rot_ke_mol_cpts = np.zeros((nmolecules, 3))  # Rotational KE of each molecule - xyz cpts
    vib_ke_mol = np.zeros(nmolecules)  # Vib KE of each molecule
    vib_ke_mol_cpts = np.zeros((nmolecules, 3))  # Vib KE of each molecule - xyz cpts
    mix1_ke_mol = np.zeros(nmolecules)  # Mixed KE of each molecule
    mix2_ke_mol = np.zeros(nmolecules)  # Mixed KE of each molecule
    mix3_ke_mol = np.zeros(nmolecules)  # Mixed KE of each molecule
    for i in range(nmolecules):  # For each molecule
        mol_id = mol_ids_unique[i]
        atom_inds = (mol_ids == mol_id)
        mol_elements = element[atom_inds]
        type_counts[i, 0] = np.sum(mol_elements == "H")
        type_counts[i, 1] = np.sum(mol_elements == "C")
        type_counts[i, 2] = np.sum(mol_elements == "O")
        natoms = mol_atom_counts[i]
        assert np.sum(atom_inds) == natoms
        msplit = mass_arr[atom_inds]
        mass_mol[i] = np.sum(msplit)
        xyzsplit = xyz_uw[atom_inds, :]
        vxyzsplit = vxyz[atom_inds, :]
        # Compute total KE of each molecule
        tot_ke_mol[i] = 0.5 * np.sum(msplit * np.sum(vxyzsplit * vxyzsplit, axis=1))
        tot_ke_mol_cpts[i] = 0.5 * np.sum(msplit[:, np.newaxis] * vxyzsplit * vxyzsplit, axis=0)
        # Compute centre of mass stats
        com[i] = np.sum(msplit[:, np.newaxis] * xyzsplit, axis=0) / mass_mol[i]
        vcom[i] = np.sum(msplit[:, np.newaxis] * vxyzsplit, axis=0) / mass_mol[i]
        # Compute angular momentum
        if natoms > 1:  # Keep ang_vel and rotational KE as 0 when only 1 atom in molecule
            angmom = np.zeros(3)
            inertia = np.zeros((3, 3))
            dxyz = xyzsplit - com[i]
            for j in range(natoms):
                angmom += msplit[j] * np.cross(dxyz[j], vxyzsplit[j])
                inertia[0, 0] += msplit[j] * (dxyz[j, 1]**2 + dxyz[j, 2]**2)
                inertia[1, 1] += msplit[j] * (dxyz[j, 0]**2 + dxyz[j, 2]**2)
                inertia[2, 2] += msplit[j] * (dxyz[j, 0]**2 + dxyz[j, 1]**2)
                inertia[0, 1] -= msplit[j] * dxyz[j, 0] * dxyz[j, 1]
                inertia[0, 2] -= msplit[j] * dxyz[j, 0] * dxyz[j, 2]
                inertia[1, 2] -= msplit[j] * dxyz[j, 1] * dxyz[j, 2]
            inertia[1, 0] = inertia[0, 1]
            inertia[2, 0] = inertia[0, 2]
            inertia[2, 1] = inertia[1, 2]
            if natoms == 2:
                rot_ke_mol[i], ang_vel[i] = rot_analysis_two_atoms(angmom, inertia, xyzsplit, check=0)
            else:
                rot_ke_mol[i], ang_vel[i] = rot_analysis(angmom, inertia, check=0)
            rot_ke_mol_cpts[i] = 0.5 * angmom * ang_vel[i]
            # Compute vib velocity
            v_vibsplit = np.zeros((natoms, 3))
            for j in range(natoms):
                v_vibsplit[j] = vxyzsplit[j] - vcom[i] - np.cross(ang_vel[i], dxyz[j])
            vib_ke_mol[i] = 0.5 * np.sum(msplit * np.sum(v_vibsplit * v_vibsplit, axis=1))
            vib_ke_mol_cpts[i] = 0.5 * np.sum(msplit[:, np.newaxis] * v_vibsplit * v_vibsplit, axis=0)
            mix1_ke_mol[i] = np.sum(msplit * np.sum(np.cross(ang_vel[i], dxyz) * vcom[i], axis=1))
            mix2_ke_mol[i] = np.sum(msplit * np.sum(vcom[i] * v_vibsplit, axis=1))
            mix3_ke_mol[i] = np.sum(msplit * np.sum(np.cross(ang_vel[i], dxyz) * v_vibsplit, axis=1))

    assert math.isclose(np.sum(mass_mol), np.sum(mass_arr))
    tr_ke_mol = 0.5 * mass_mol * np.sum(vcom * vcom, axis=1)  # Translational KE of each molecule
    tr_ke_mol_cpts = 0.5 * mass_mol[:, np.newaxis] * vcom * vcom  # Translational KE of each molecule - xyz cpts
    # Test xyz cpts of KE sum to total
    np.testing.assert_allclose(np.sum(tot_ke_mol_cpts, axis=1), tot_ke_mol, atol=1e-6)
    np.testing.assert_allclose(np.sum(tr_ke_mol_cpts, axis=1), tr_ke_mol, atol=1e-6)
    np.testing.assert_allclose(np.sum(rot_ke_mol_cpts, axis=1), rot_ke_mol, atol=1e-6)
    np.testing.assert_allclose(np.sum(vib_ke_mol_cpts, axis=1), vib_ke_mol, atol=1e-6)
    # Convert energies to kcal/mol
    tot_ke_mol *= CONV_TO_KCAL_MOL
    tot_ke_mol_cpts *= CONV_TO_KCAL_MOL
    tr_ke_mol *= CONV_TO_KCAL_MOL
    tr_ke_mol_cpts *= CONV_TO_KCAL_MOL
    rot_ke_mol *= CONV_TO_KCAL_MOL
    rot_ke_mol_cpts *= CONV_TO_KCAL_MOL
    vib_ke_mol *= CONV_TO_KCAL_MOL
    vib_ke_mol_cpts *= CONV_TO_KCAL_MOL
    mix1_ke_mol *= CONV_TO_KCAL_MOL
    mix2_ke_mol *= CONV_TO_KCAL_MOL
    mix3_ke_mol *= CONV_TO_KCAL_MOL
    prop = (tr_ke_mol + rot_ke_mol + vib_ke_mol + mix1_ke_mol + mix2_ke_mol + mix3_ke_mol) / tot_ke_mol
    assert np.all(np.isclose(prop, 1.0))
    # Analyze size of mixed KE terms
    tol_size = 1e-6
    if np.any(mix1_ke_mol / tot_ke_mol > tol_size):
        print(f"{mix1_ke_mol = }")
    if np.any(mix2_ke_mol / tot_ke_mol > tol_size):
        print(f"{mix2_ke_mol = }")
    if np.any(mix3_ke_mol / tot_ke_mol > tol_size):
        print(f"{mix3_ke_mol = }")
    counts_type = np.column_stack((mol_atom_counts, type_counts))
    return tot_ke_mol_cpts, tr_ke_mol_cpts, rot_ke_mol_cpts, vib_ke_mol_cpts, counts_type

def analyze(dump_file, save=0):
    """ Analyze kinetic energies """
    ntimesteps, natoms, box_len = overview(dump_file)
    mass_arr = mass_array(dump_file, natoms)
    # Initialise per-atom arrays
    atom_ids = np.zeros(natoms, dtype=int)
    element = np.zeros(natoms, dtype=object)
    q = np.zeros(natoms)
    mol_ids = np.zeros(natoms, dtype=int)
    xyz = np.zeros((natoms, 3))
    vxyz = np.zeros((natoms, 3))
    ixyz = np.zeros((natoms, 3), dtype=int)
    # Initialise arrays over time - slice at end to account for every (when skipping timesteps)
    timesteps = np.zeros(ntimesteps, dtype=int)
    all_ke = np.zeros((ntimesteps, 12))  # Kinetic energy of all species - total, trans, rot, vib
    fuel_ke = np.zeros((ntimesteps, 12))  # Kinetic energy of C12H26 - total, trans, rot, vib
    nfuel = np.zeros(ntimesteps, dtype=int)
    oxy_ke = np.zeros((ntimesteps, 12))  # Kinetic energy of oxygen molecules - total, trans, rot, vib
    noxy = np.zeros(ntimesteps, dtype=int)
    dof_counts = np.zeros((ntimesteps, 3), dtype=int)  # Trans, rot, vib dofs of all species
    save_count = 0
    every = 10000
    # Loop through timesteps
    f = open(dump_file, "r")
    for nt in range(ntimesteps):
        # Domain data
        f.readline()
        timestep = int(f.readline())
        for _ in range(7):
            f.readline()
        if timestep % every == 0:
            if timestep % 10000 == 0:
                print(f"{timestep = }")
            # Read in atom data
            for i in range(natoms):
                lsplit = f.readline().split()
                atom_ids[i] = int(lsplit[0])
                element[i] = lsplit[1] 
                q[i] = float(lsplit[2])
                mol_ids[i] = int(lsplit[3])
                xyz[i, 0] = float(lsplit[4])
                xyz[i, 1] = float(lsplit[5])
                xyz[i, 2] = float(lsplit[6])
                vxyz[i, 0] = float(lsplit[7])
                vxyz[i, 1] = float(lsplit[8])
                vxyz[i, 2] = float(lsplit[9])
                ixyz[i, 0] = int(lsplit[10])
                ixyz[i, 1] = int(lsplit[11])
                ixyz[i, 2] = int(lsplit[12])
            # Check atom ids
            np.testing.assert_equal(atom_ids, np.linspace(1, natoms, natoms, dtype=int))
            # Analyze molecules
            xyz_uw = xyz + box_len * ixyz  # Unwrapped coords
            tot_ke_mol_cpts, tr_ke_mol_cpts, rot_ke_mol_cpts, vib_ke_mol_cpts, counts_type \
                    = mol_analysis(mol_ids, element, mass_arr, xyz_uw, vxyz)
            # Check total energies match
            total_ke = 0.5 * np.sum(mass_arr[:, np.newaxis] * vxyz * vxyz, axis=0) * CONV_TO_KCAL_MOL
            np.testing.assert_allclose(np.sum(tot_ke_mol_cpts, axis=0), total_ke, atol=1e-6)
            # Analysis for saving
            timesteps[save_count] = timestep
            # Total
            all_ke[save_count, 0:3] = np.sum(tot_ke_mol_cpts, axis=0)
            all_ke[save_count, 3:6] = np.sum(tr_ke_mol_cpts, axis=0)
            all_ke[save_count, 6:9] = np.sum(rot_ke_mol_cpts, axis=0)
            all_ke[save_count, 9:12] = np.sum(vib_ke_mol_cpts, axis=0)
            # C12H26 stats
            tmp = np.logical_and(counts_type[:, 1] == 26, counts_type[:, 2] == 12)  # 26 H and 12 C
            hc_loc = np.logical_and(tmp, counts_type[:, 0] == 38)  # Total atoms = 38
            fuel_ke[save_count, 0:3] = np.sum(tot_ke_mol_cpts[hc_loc], axis=0)
            fuel_ke[save_count, 3:6] = np.sum(tr_ke_mol_cpts[hc_loc], axis=0)
            fuel_ke[save_count, 6:9] = np.sum(rot_ke_mol_cpts[hc_loc], axis=0)
            fuel_ke[save_count, 9:12] = np.sum(vib_ke_mol_cpts[hc_loc], axis=0)
            nfuel[save_count] = np.sum(hc_loc)
            # O2 stats
            oxy_loc = np.logical_and(counts_type[:, 0] == 2, counts_type[:, 3] == 2)  # 2 O and 2 atoms
            oxy_ke[save_count, 0:3] = np.sum(tot_ke_mol_cpts[oxy_loc], axis=0)
            oxy_ke[save_count, 3:6] = np.sum(tr_ke_mol_cpts[oxy_loc], axis=0)
            oxy_ke[save_count, 6:9] = np.sum(rot_ke_mol_cpts[oxy_loc], axis=0)
            oxy_ke[save_count, 9:12] = np.sum(vib_ke_mol_cpts[oxy_loc], axis=0)
            noxy[save_count] = np.sum(oxy_loc)
            # All dof counts
            nmolecules = len(tot_ke_mol_cpts)
            monoatomic = np.sum(counts_type[:, 0] == 1)
            diatomic = np.sum(counts_type[:, 0] == 2)
            # Consider linear molecules to be diatomic and carbon dioxide
            co2_loc = np.logical_and(np.logical_and(counts_type[:, 2] == 1, counts_type[:, 3] == 2),
                                     counts_type[:, 0] == 3)
            co2_count = np.sum(co2_loc)
            # Counts of linear and nonlinear molecules
            linear = diatomic + co2_count
            nonlinear = nmolecules - monoatomic - linear
            # Number of atoms in nonlinear molecules
            natoms_linear = 2 * diatomic + 3 * co2_count
            natoms_nonlinear = np.sum(counts_type[counts_type[:, 0] > 2, 0]) - 3 * co2_count
            dof_counts[save_count, 0] = 3 * nmolecules  # Translation
            dof_counts[save_count, 1] = 2 * linear + 3 * nonlinear  # Rotation
            dof_counts[save_count, 2] = 3 * (natoms_linear + natoms_nonlinear) \
                                        - 5 * linear - 6 * nonlinear  # Vibration
            assert np.sum(dof_counts[save_count]) == 3 * natoms, f"{np.sum(dof_counts) = }"

            save_count += 1
        else:
            for _ in range(natoms):
                f.readline()

    f.close()

    if save:
        # Slice np arrays to account for timesteps skipped
        timesteps = timesteps[:save_count]
        all_ke = all_ke[:save_count, :]
        fuel_ke = fuel_ke[:save_count, :]
        nfuel = nfuel[:save_count]
        oxy_ke = oxy_ke[:save_count, :]
        noxy = noxy[:save_count]
        dof_counts = dof_counts[:save_count, :]
        ke_dir = os.path.join(os.path.dirname(dump_file), "ke_analysis")
        if not os.path.exists(ke_dir):
            os.makedirs(ke_dir)
        def save_ke(ke, n, name, header_extra):
            """ Save ke stats to file """
            header_start = "Timestep, total KE (x y z), trans KE (x y z), rot KE (x y z), vib KE (x y z), "
            header = header_start + header_extra + ". All KE in kcal/mol"
            fmt = "%8.u %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f"
            save_arr = np.column_stack((timesteps, ke, n))
            if np.shape(save_arr)[1] == 14:
                fmt += " %6.u"
            elif np.shape(save_arr)[1] == 16:
                fmt += " %6.u %6.u %6.u"
            else:
                raise Exception(f"Unexpected {np.shape(save_arr) = }")
            save_file = os.path.join(ke_dir, f"{name}.txt")
            np.savetxt(save_file, save_arr, header=header, fmt=fmt)
            print(f"Saved {save_file = }")

        save_ke(all_ke, dof_counts, "all_stats", header_extra="dof counts (tr, rot, vib)")
        save_ke(fuel_ke, nfuel, "fuel_stats", header_extra="no. of fuel molecules")
        save_ke(oxy_ke, noxy, "oxy_stats", header_extra="no. of oxygen molecules")

if __name__ == "__main__":
    CONV_TO_KCAL_MOL = 1e4 / 4.184

    analyze("../pmh_combustion/qeq_run_nhoover_temp2000/ef_0.0/run_1/dump.run", save=1)
