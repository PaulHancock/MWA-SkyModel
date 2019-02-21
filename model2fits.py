#! /usr/bin/env python

from astropy.table import Table
from astropy.coordinates import SkyCoord
import os

__author__ = 'Paul Hancock'
__date__ = '2019-02-21'


class Base(object):
    def __init__(self, string):
        self.value = string

    def __repr__(self):
        return self.value


class Name(Base):
    def __init__(self, string):
        Base.__init__(self, string.replace('"', '').replace("'", ""))


class SrcType(Base):
    def __init__(self, string):
        if str.lower(string) == 'gaussian':
            self.value = 'gaussian'
        else:
            self.value = 'point'


class Position(Base):
    def __init__(self, string):
        Base.__init__(self, string)
        self.ra_str, self.dec_str = string.split()
        c = SkyCoord(self.ra_str, self.dec_str)
        self.ra = c.ra.degree
        self.dec = c.dec.degree

    def totable(self):
        c = SkyCoord(self.ra_str, self.dec_str)
        # convert to hmsdms style
        s = c.to_string(style='hmsdms')
        # use colon separators
        s = s.replace('h', ':').replace('m', ':').replace('s', '').replace('d', ':')
        return s


class Shape(Base):
    def __init__(self, string):
        Base.__init__(self, string)
        self.a, self.b, self.pa = list(map(float, string.split()))


class Spec(Base):
    def __init__(self, alpha, beta=0):
        Base.__init__(self, "{0} {1}".format(alpha, beta))
        self.alpha, self.beta = alpha, beta

    def __repr__(self):
        r = 'spectral-index {\n'
        r += '  alpha {0}\n'.format(self.alpha)
        r += '  beta {0}\n'.format(self.beta)
        r += '}'
        return r


class Measurement(object):
    def __init__(self, freq, flux):
        self.freq, self.freq_units = freq.split(' ')
        self.freq = float(self.freq)
        self.flux_units = flux.split(' ')[0]
        self.I, self.Q, self.U, self.V = list(map(float, flux.split(' ')[1:]))

    def __repr__(self):
        r = 'measurement {\n'
        r += '  frequency {0} {1}\n'.format(self.freq, self.freq_units)
        r += '  fluxdensity {0} {1} {2} {3} {4}\n'.format(self.flux_units, self.I, self.Q, self.U, self.V)
        r += '}'
        return r


class Component(object):
    def __init__(self, children):
        self._type = None
        self.shape = None
        self.position = None
        self.spec = Spec(0,0)
        self.measurement = []
        c = iter(children)
        while True:
            try:
                l = c.next().strip()
                if len(l) < 1 or l.startswith('}'):
                    break
                print(l)
                key, val = l.split(' ', 1)
                print("processing {0} = {1}".format(key, val))
                if key.startswith('type'):
                    self._type = SrcType(val)
                elif key.startswith('shape'):
                    self.shape = Shape(val)
                elif key.startswith('position'):
                    self.position = Position(val)
                elif key.startswith('measurement'):
                    _, freq = c.next().strip().split(' ', 1)
                    _, flux = c.next().strip().split(' ', 1)
                    print("making measurement with |{0}| / |{1}|".format(freq, flux))
                    self.measurement.append(Measurement(freq, flux))
                    c.next()  # skip the closing brace
                elif key.startswith('spectral'):
                    _, alpha = c.next().strip().split(' ', 1)
                    _, beta = c.next().strip().split(' ', 1)
                    print("making spectral-index with |{0}| / |{1}|".format(alpha, beta))
                    self.spec = Spec(alpha, beta)
                    c.next()  # skip the closing brace

                else:
                    print("{0} not recognized".format(key))
            except StopIteration:
                break

    def __repr__(self):
        r = 'component {\n'
        r += '  type {0}\n'.format(self._type)
        if self._type == 'gaussian':
            r += '  shape {0}\n'.format(self.shape)
        r += '  position {0}\n'.format(self.position)
        if self.spec is not None:
            for i in self.spec.__repr__().split('\n'):
                r += '  {0}\n'.format(i)
        for m in self.measurement:
            for i in m.__repr__().split('\n'):
                r += '  {0}\n'.format(i)
        r += '}'
        return r


class Source(object):
    def __init__(self, data):
        self.name = None
        self.components = []
        c = iter(data)
        while True:
            try:
                l = c.next().strip()
                if len(l) < 1 or l.startswith('}'):
                    break
                print(l)
                key, val = l.split(' ', 1)
                print("processing {0} = {1}".format(key, val))
                if key.startswith('name'):
                    self.name = Name(val)
                elif key.startswith('component'):
                    # find the closing brace
                    brace = 1
                    children = []
                    while brace > 0:
                        l = c.next().strip()
                        if l.startswith('}'):
                            brace -= 1
                        elif l[-1] == '{':
                            brace += 1
                        children.append(l)
                    print("making component")
                    self.components.append(Component(children))
                else:
                    print("{0} not recognized".format(key))
            except StopIteration:
                break

    def __repr__(self):
        r = 'source {\n'
        r += '  name "{0}"\n'.format(self.name)
        for c in self.components:
            for i in c.__repr__().split('\n'):
                r += '  {0}\n'.format(i)
        r += '}'
        return r

    def as_table(self):
        basename = self.name
        rows = []
        for c in self.components:
            # name, ra_str, dec_str, a_wide, b_wide, pa_wide, peak_flux_wide
            ra_str, dec_str = c.position.totable().split(' ', 1)
            if c.shape is None:
                a, b, pa = 0., 0., 0.
            else:
                a, b, pa = c.shape.a, c.shape.b, c.shape.pa
            row = [basename,
                   c.position.ra, c.position.dec, ra_str, dec_str,
                   a, b, pa,
                   c.measurement[0].freq, c.measurement[0].I,
                   c.spec.alpha]
            rows.append(row)
        tab = Table(rows=rows,
                    names=('Name', 'ra', 'dec', 'ra_str', 'dec_str', 'a', 'b', 'pa', 'freq', 'peak_flux', 'alpha'),
                    dtype=(str, float, float, str, str, float, float, float, float, float, float))
        return tab


if __name__ == '__main__':
    f = 'CenA_core_wsclean_model.txt'
    lines = open(f).readlines()
    src = Source(lines)
    out = 'Skymodel.fits'
    if os.path.exists(out):
        os.remove(out)
    src.as_table().write(out)
    pass
