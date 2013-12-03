import sys

if len(sys.argv) < 2:
    print 'Usage:', sys.argv[0], 'element.txt'
    exit(-1)

f = open(sys.argv[1])
lines = f.readlines()

elements = {}

current = 1
for line in lines:
    # slip comments
    if line[0] == '#':
        continue

    columns = ' '.join(line[:-1].split()).split(' ')
    if len(columns) != 9:
        print 'Error: incorrect number of columns for line:'
        print line
        exit(-1)

    elements[current] = columns

    current += 1

    
#print elements

print '#ifndef HELIUM_ELEMENT_H'
print '#define HELIUM_ELEMENT_H'
print ''
print '#include <string>'
print ''
print 'namespace Helium {'
print ''
print '  /**'
print '   * @brief Element properties.'
print '   */'
print '  class Element'
print '  {'
print '    public:'

print '      /**'
print '       * @brief Get the element symbol for an atom.'
print '       *'
print '       * @param element The atom element number.'
print '       *'
print '       * @return The element symbol.'
print '       */'
print '      static std::string symbol(int element)'
print '      {'
print '        switch (element) {'
for element in elements.keys():
    print '          case ' + str(element) + ':'
    print '            return "' + str(elements[element][0]) + '";'
print '          default:'
print '            return "Xx";'
print '        }'
print '      }'

print ''

print '      /**'
print '       * @brief Get the valence for an atom.'
print '       *'
print '       * @param element The atom element number.'
print '       *'
print '       * @return The average mass.'
print '       */'
print '      static int averageMass(int element)'
print '      {'
print '        switch (element) {'
for element in elements.keys():
    print '          case ' + str(element) + ':'
    print '            return ' + str(elements[element][1]) + ';'
print '          default:'
print '            return 0;'
print '        }'
print '      }'

print ''

print '      /**'
print '       * @brief Check if hydrogens should be added for this atom.'
print '       *'
print '       * @param element The atom element number.'
print '       *'
print '       * @return True if hydrogens should be added.'
print '       */'
print '      static bool addHydrogens(int element)'
print '      {'
print '        switch (element) {'
for element in elements.keys():
    if elements[element][3] == 'Yes':
        print '          case ' + str(element) + ':'
print '            return true;'

print '          default:'
print '            return false;'
print '        }'
print '      }'

print ''

print '      /**'
print '       * @brief Get the valence for an atom'
print '       *'
print '       * @param element The atom element number.'
print '       * @param charge The atom charge.'
print '       * @param degree The current atom degree.'
print '       *'
print '       * @return The valence.'
print '       */'
print '      static int valence(int element, int charge, int degree)'
print '      {'
print '        if (charge < -2 || charge > 2)'
print '          return 0;'
print '        switch (element) {'
for element in elements.keys():
    print '          case ' + str(element) + ':'
    allDashes = True
    for charge in range(5):
        if elements[element][4 + charge] != '-':
            allDashes = False
    if allDashes:
        print '            return 0;'
        continue
    print '            switch (charge) {'
    for charge in range(5):
        print '              case ' + str(charge - 2) + ':'
        if elements[element][4 + charge] == '-':
            print '                return 0;'
        else:
            valences = elements[element][4 + charge].split(',')
            for valence in valences:
                if len(valence) == 2:
                    continue
                print '                if (degree < ' + valence + ')'
                print '                  return ' + valence + ';'
            print '                return 0;'
    print '            }'
print '          default:'
print '            return 0;'
print '        }'
print '      }'


print '  };'
print ''
print '}'
print ''
print '#endif'
