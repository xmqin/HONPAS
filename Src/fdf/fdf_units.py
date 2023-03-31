#
# Small script to auto-generate the fdf_convfac
# parameters used for conversion between units.
# This is merely due to the easy conversion of data
# in python, relative to manually adding elements in
# FORTRAN
from __future__ import print_function

from collections import OrderedDict

units = OrderedDict()

units['mass'] = OrderedDict()
mass = units['mass']
mass['g'] = '1.d-3'
mass['kg'] = '1.d0'
mass['amu'] = '1.66054d-27'

units['length'] = OrderedDict()
length = units['length']
length['m'] = '1.d0'
length['cm'] = '1.d-2'
length['nm'] = '1.d-9'
length['pm'] = '1.d-12'
length['Ang'] = '1.d-10'
length['Bohr'] = '0.529177d-10'

units['energy'] = OrderedDict()
energy = units['energy']
energy['J'] = '1.d0'
energy['kJ'] = '1.d3'
energy['erg'] = '1.d-7'
energy['meV'] = '1.60219d-22'
energy['eV'] = '1.60219d-19'
energy['mRy'] = '2.17991d-21'
energy['Ry'] = '2.17991d-18'
energy['mHa'] = '4.35982d-21'
energy['mHartree'] = '4.35982d-21'
energy['Ha'] = '4.35982d-18'
energy['Hartree'] = '4.35982d-18'
energy['K'] = '1.38066d-23'
energy['Kelvin'] = '1.38066d-23'
energy['kcal/mol'] = '6.94780d-21'
energy['kJ/mol'] = '1.6606d-21'
energy['Hz'] = '6.6262d-34'
energy['THz'] = '6.6262d-22'
energy['cm-1'] = '1.986d-23'
energy['cm^-1'] = '1.986d-23'
energy['cm**-1'] = '1.986d-23'

units['time'] = OrderedDict()
time = units['time']
time['s'] = '1.d0'
time['ns'] = '1.d-9'
time['ps'] = '1.d-12'
time['fs'] = '1.d-15'
time['min'] = '60.d0'
time['mins'] = '60.d0'
time['hour'] = '3600.d0'
time['hours'] = '3600.d0'
time['day'] = '86400.d0'
time['days'] = '86400.d0'

units['force'] = OrderedDict()
force = units['force']
force['N'] = '1.d0'
force['eV/Ang'] = '1.60219d-9'
force['Ry/Bohr'] = '4.11943d-8'
force['Ha/Bohr'] = '2.059715d-08'

units['pressure'] = OrderedDict()
pressure = units['pressure']
pressure['Pa'] = '1.d0'
pressure['GPa'] = '1.d9'
pressure['atm'] = '1.01325d5'
pressure['bar'] = '1.d5'
pressure['Mbar'] = '1.d11'
pressure['eV/Ang**3'] = '1.60219d11'
pressure['eV/Ang^3'] = '1.60219d11'
pressure['Ry/Bohr**3'] = '1.47108d13'
pressure['Ry/Bohr^3'] = '1.47108d13'

units['charge'] = OrderedDict()
charge = units['charge']
charge['c'] = '1.d0'
charge['e'] = '1.602177d-19'

units['dipole'] = OrderedDict()
dipole = units['dipole']
dipole['c*m'] = '1.d0'
dipole['D'] = '3.33564d-30'
dipole['Debye'] = '3.33564d-30'
dipole['e*Bohr'] = '8.47835d-30'
dipole['e*Ang'] = '1.602177d-29'

units['mominert'] = OrderedDict()
mominert = units['mominert']
mominert['kg*m**2'] = '1.d0'
mominert['Ry*fs**2'] = '2.17991d-48'

units['efield'] = OrderedDict()
efield = units['efield']
efield['V/m'] = '1.d0'
efield['V/nm'] = '1.d9'
efield['V/Ang'] = '1.d10'
efield['V/Bohr'] = '1.8897268d10'
efield['Ry/Bohr/e'] = '2.5711273d11'
efield['Ha/Bohr/e'] = '5.1422546d11'
efield['Har/Bohr/e'] = '5.1422546d11'

units['angle'] = OrderedDict()
angle = units['angle']
angle['deg'] = '1.d0'
angle['rad'] = '5.72957795d1'

units['torque'] = OrderedDict()
torque = units['torque']
torque['meV/deg'] = '1.0d-3'
torque['meV/rad'] = '1.745533d-5'
torque['eV/deg'] = '1.0d0'
torque['eV/rad'] = '1.745533d-2'
torque['mRy/deg'] = '13.6058d-3'
torque['mRy/rad'] = '0.237466d-3'
torque['Ry/deg'] = '13.6058d0'
torque['Ry/rad'] = '0.237466d0'


def check_ambiguity(units):
    for field1 in units:
        unit1 = units[field1]
        for name1 in unit1:
            i = 0
            lst = []
            for field2 in units:
                unit2 = units[field2]
                for name2 in unit2:
                    if name1.lower() == name2.lower():
                        i += 1
                        if i > 1:
                            lst.append((field1, name1, field2, name2))
            if i != 1:
                print(lst)
                raise ValueError('Ambiguity in units, several has the same name.')

check_ambiguity(units)

def max_length(units):
    fl = 1
    ul = 1
    for field in units:
        fl = max(fl, len(field))
        for unit in units[field]:
            ul = max(ul, len(unit))
    return fl, ul

def length(units):
    l = 0
    for field in units:
        l = l + len(units[field])
    return l


dimm_l, name_l = max_length(units)
max_units = 10

# Now write everything
ind = ' ' * 6
print(ind + 'integer(ip), parameter :: nu = {}'.format(length(units)))
print(ind + 'character({}) :: dimm(nu)'.format(dimm_l))
print(ind + 'character({}) :: name(nu)'.format(name_l))
print(ind + 'real(dp) :: unit(nu)')

fmt_s = "'{0:" + str(dimm_l) + "s}', '{1:" + str(name_l) + "s}', {2:s}"
# Double precision has (up to) 17 significant digits.
fmt_f = "'{0:" + str(dimm_l) + "s}', '{1:" + str(name_l) + "s}', {2:<.17e}_dp"

def get_line(field, unit, val):
    if isinstance(val, float):
        return fmt_f.format(field, unit, val)
    else:
        return fmt_s.format(field, unit, val)

N = 0
for field in units:
    # Number of units for this field
    nunits = len(units[field])
    
    iind = ind + ' ' * 4
    for i, unit in enumerate(units[field]):
        if i % max_units == 0:
            if i > 0:
                N += n
            n = min(len(units[field])-i, max_units)
            print(ind + 'data (dimm(iu), name(iu), unit(iu), iu={}, {}) / &'.format(N+1, N+n))
        if i % max_units == n - 1:
            # End with '/
            print(iind + get_line(field, unit.lower(), units[field][unit]) + ' /')
        else:
            print(iind + get_line(field, unit.lower(), units[field][unit]) + ', &')
    print()
    N += n
