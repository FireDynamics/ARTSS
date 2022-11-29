import sys
import numpy as np
import matplotlib.pyplot as plt
from ARTSS import XML, Domain


def plot(x, y, z, show=False):
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(x, y, z)
    plt.xticks(np.unique(x))
    plt.xlabel('x')
    plt.yticks(np.unique(y))
    plt.ylabel('y')
    ax.set_zticks(np.unique(z))
    ax.set_zlabel('z')
    if show:
        plt.show()
    else:
        plt.savefig('stern.pdf')
        plt.close()


def main(file_name: str):
    xml = XML(file_name)
    xml.read_xml()
    domain = Domain(domain_param=xml.domain, obstacles=xml.obstacles, enable_computational_domain=xml.computational_domain)
    domain.print_info()
    #    domain.print_debug()
    while True:
        text = input('index, coordinates "i j k", combo: ')
        if text.startswith('combo'):
            coords = input('i1 i2 j1 j2 k1 k2: ')
            array = coords.split(' ')
            indices = [int(a) for a in array]
            x = []
            y = []
            z = []
            for i in indices[:2]:
                for j in indices[2:4]:
                    for k in indices[4:]:
                        index = domain.calculate_index(i, j, k)
                        print(domain.get_type(index))
                        x.append(i)
                        y.append(j)
                        z.append(k)
            answer = input('draw? [y/N] ')
            answer.lower()
            if answer in ('y', 'yes'):
                data = input('original point: ')
                data = data.split(' ')
                if len(data) > 2:
                    x.append(int(data[0]))
                    y.append(int(data[1]))
                    z.append(int(data[2]))
                plot(x, y, z, show=True)
        elif text.startswith('neighbor'):
            index_original = domain.calculate_index(int(array[0]), int(array[1]), int(array[2]))
            print('orig.\t' + domain.get_type(int(index_original)))
            index_front = domain.calculate_index(int(array[0]), int(array[1]), int(array[2]) - 1)
            print('front\t' + domain.get_type(int(index_front)))
            index_back = domain.calculate_index(int(array[0]), int(array[1]), int(array[2]) + 1)
            print('back\t' + domain.get_type(int(index_back)))
            index_bottom = domain.calculate_index(int(array[0]), int(array[1]) - 1, int(array[2]))
            print('bottom\t' + domain.get_type(int(index_bottom)))
            index_top = domain.calculate_index(int(array[0]), int(array[1]) + 1, int(array[2]))
            print('top\t' + domain.get_type(int(index_top)))
            index_left = domain.calculate_index(int(array[0]) - 1, int(array[1]), int(array[2]))
            print('left\t' + domain.get_type(int(index_left)))
            index_right = domain.calculate_index(int(array[0]) + 1, int(array[1]), int(array[2]))
            print('right\t' + domain.get_type(int(index_right)))
        else:
            array = text.split(' ')
            if len(array) == 3:
                index = domain.calculate_index(int(array[0]), int(array[1]), int(array[2]))
                print(domain.get_type(index))
            else:
                print(domain.get_type(int(text)))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('filename is mandatory')
    main(sys.argv[1])
