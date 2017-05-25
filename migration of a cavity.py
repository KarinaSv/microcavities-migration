#Программа, для описания миграции микрополости в полимерном композиционном материале

import math
import random

#Константы
r = 20 #радиус микрополости
R = 50 #радиус наполнителя
mol_R = 11 #радиус макромолекул полимера
border = 500 #граница области
Try = 100  # количество попыток движений микропоры
volume_cavity = 4 * math.pi * r ** 3 / 3 #объем полости
distance = 200 #расстояние между наполнителями
count = 10 ** 3

class Sphere:
    def __init__(self, x, y, z, radius):
        self.x = x
        self.y = y
        self.z = z
        self.radius = radius
        self.volume = 4 * math.pi * self.radius ** 3 / 3

#Розыгрыш начальных координат микропоры
def micropore(border, microradius):
    x = (border - 2 * microradius) * random.random()
    y = (border - 2 * microradius) * random.random()
    z = (border - 2 * microradius) * random.random()
    return x, y, z

#Розыгрыш координат конца вектора свободного пробега
def line(border, R):
    x = border * random.random()
    y = border * random.random()
    z = (border + 2 * R) * random.random()
    return x, y, z

#Распределение наполнителей в заданном объеме
def get_fillers(distance, border, filler_radius):
    filler_list = []
    x = filler_radius
    while True:
        if (x < border):
            y = filler_radius
            while True:
                if (y < border):
                    z = filler_radius
                    while True:
                        if (z < border):
                            new_sphere = Sphere(x, y, z, filler_radius)
                            filler_list.append(new_sphere)
                            z += distance
                        else:
                            break
                    y += distance
                else:
                    break
            x += distance
        else:
            break

    return filler_list

#Розыгрыш координат макромолекул полимера
def molecule(border, mol_radius):
    x = (border - 2 * mol_radius) * random.random() + mol_radius
    y = (border - 2 * mol_radius) * random.random() + mol_radius
    z = (border - 2 * mol_radius) * random.random() + mol_radius
    return x, y, z

#Распределение макромолекул полимера в заднном объеме
def get_molecules(count, border, mol_radius, filler_list):
    mol_list = []
    max_number_of_tries = 10 ** 5
    number_of_tries = 0
    while len(mol_list) < count:
        x, y, z = molecule(border, mol_radius)
        new_sphere = Sphere(x, y, z, mol_radius)
        has_crossed = False

        for fil in mol_list:
            if abs(fil.x - new_sphere.x) < 2 * mol_radius or \
               abs(fil.y - new_sphere.y) < 2 * mol_radius or \
               abs(fil.z - new_sphere.z) < 2 * mol_radius:
                has_crossed = True
                break

        for filler in filler_list:
            if abs(filler.x - new_sphere.x) < (filler.radius + new_sphere.radius) and abs(filler.y - new_sphere.y) < (
                filler.radius + new_sphere.radius) and abs(filler.z - new_sphere.z) < (filler.radius + new_sphere.radius):
                has_crossed = True
                break

        if has_crossed:
            number_of_tries += 1
            if number_of_tries >= max_number_of_tries:
                break
        else:
            mol_list.append(new_sphere)
            number_of_tries = 0

    return mol_list

#Проверка столкновения траектории полости с элементами материала
def is_collision(x0, y0, z0, x, y, z, xk, yk, zk, filler_radius):
    a = (xk - x0) ** 2 + (yk - y0) ** 2 + (zk - z0) ** 2
    b = 2 * ((xk - x0) * (x0 - x) + (yk - y0) * (y0 - y) + (zk - z0) * (z0 - z))
    с = x ** 2 + y ** 2 + z ** 2 + x0 ** 2 + y0 ** 2 + z0 **2 - 2 * (x * x0 + y * y0 + z * z0) - filler_radius ** 2
    D = b ** 2 - 4 * a * с
    return D, a, b, с

#Обработка столкновения
def collision(a, b, c, x0, y0, z0, xk, yk, zk):
    t1 = (-b + math.sqrt(b ** 2 - 4 * a *c)) / (2 * a)
    t2 = (-b - math.sqrt(b ** 2 - 4 * a *c)) / (2 * a)

    #Расчет точек пересечения прямой и сферы
    def cross(start, end, t):
        new_end = start + (end - start) * t
        return new_end

    xk1 = cross(x0, xk, t1) #расчет первой точки пересения
    yk1 = cross(y0, yk, t1)
    zk1 = cross(z0, zk, t1)

    xk2 = cross(x0, xk, t2) #расчет второй точки пересечения
    yk2 = cross(y0, yk, t2)
    zk2 = cross(z0, zk, t2)

    #Проверка первого пересечения со сферой
    if abs(xk1 - x0) < abs(xk2 - x0):
        xk = xk1
        yk = yk1
        zk = zk1
    else:
        xk = xk2
        yk = yk2
        zk = zk2

    return xk, yk, zk

#Расчет направляющих косинусов
def cos(x0, y0, z0, xk, yk, zk):
    abs_l = math.sqrt((xk - x0) ** 2 + (yk - y0) ** 2 + (zk - z0) ** 2)
    cos_alfa = (xk - x0) / abs_l
    cos_betta = (yk - y0) / abs_l
    cos_gamma = (zk - z0) / abs_l
    return cos_alfa, cos_betta, cos_gamma

#Расчет координат центра микрополости
def new_koord(xk, yk, zk, H, cos_alfa, cos_betta, cos_gamma):
    def new_centre(H, cos, end):
        centre = end - H * cos
        return centre

    x0 = new_centre(H, cos_alfa, xk)
    y0 = new_centre(H, cos_betta, yk)
    z0 = new_centre(H, cos_gamma, zk)

    return x0, y0, z0

#Расчет расстояния от точки пересечения до нового центра полости
def S_for_new_koord(x_cross, y_cross, z_cross, x0, y0, z0, x_sphere, y_sphere, z_sphere, sphere_radius, microradius):
    S_centre_to_cross = math.sqrt((x_cross - x0)**2 + (y_cross - y0)**2 + (z_cross - z0)**2)
    S_double_centre = math.sqrt((x_sphere - x0)**2 + (y_sphere - y0)**2 + (z_sphere - z0)**2)
    angle_centre_cross_centre = math.acos((sphere_radius ** 2 + S_centre_to_cross ** 2 - S_double_centre ** 2) / (2 * sphere_radius * S_centre_to_cross))
    angle_new_centre = math.asin((sphere_radius * math.sin(angle_centre_cross_centre) / (sphere_radius + microradius)))
    angle = math.pi - angle_new_centre - angle_centre_cross_centre
    H = (sphere_radius + microradius) * math.sin(angle) / math.sin(angle_centre_cross_centre)
    return (H)

#Присвоение старых координат полости другим значениям (для подсчета пути)
def old_koord(x0, y0, z0):
    old_x0 = x0
    old_y0 = y0
    old_z0 = z0
    return old_x0, old_y0, old_z0

def main():
    file_micropore = open('Micropore.txt', 'w')
    file_filler = open('Filler.txt', 'w')
    file_arrow = open('Arrow.txt', 'w')
    file_capture_micropore = open('CaptureMicropore.txt', 'w')
    file_molecules = open('Molecules.txt', 'w')

    #Cоздание наполнителей
    filler_list = get_fillers(distance, border, R)
    for fil in filler_list:
        file_filler.write(str(fil.x) + '   ' + str(fil.y) + '   ' + str(fil.z) + '\n')  # запись в файл координат центра частицы наполнителя

    #Cоздание макромолекул
    mol_list = get_molecules(count, border, mol_R, filler_list)

    #Запись координат макромолекул
    for mol in mol_list:
        file_molecules.write(str(mol.x) + '   ' + str(mol.y) + '   ' + str(mol.z) + '\n')

    #Проверка на пересечение полости с наполнителями и макромолекулами изначально
    intersection = True
    x0, y0, z0 = 0, 0, 0
    while intersection:
        intersection = False
        x0, y0, z0 = micropore(border, r)
        for fil in filler_list:
            if abs(fil.x - x0) < (R + r) and abs(fil.y - y0) < (R + r) and abs(fil.z - z0) < (R + r):
                intersection = True
                break
        for sphere in mol_list:
            if abs(sphere.x - x0) < (mol_R + r) and abs(sphere.y - y0) < (mol_R + r) and abs(sphere.z - z0) < (mol_R + r):
                intersection = True
                break
    file_micropore.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')

    Way = 0  #Путь полости
    cross_number_filler = 1000
    cross_number_mol = 1000
    capture = 0
    for i in range(1, Try):
        xk, yk, zk = line(border, r)
        for index, fil in enumerate(filler_list):
            cross_filler = 10 ** 4
            if index != cross_number_filler:
                D, a, b, c = is_collision(x0, y0, z0, fil.x, fil.y, fil.z, xk, yk, zk, R)
                if D >= 0:  # если есть пересечение
                    xk_f, yk_f, zk_f = collision(a, b, c, x0, y0, z0, xk, yk, zk)
                    cross_filler = math.sqrt((x0 - xk_f) ** 2 + (y0 - yk_f) ** 2 + (z0 - zk_f) ** 2)  #Расстояние от полости до точки пересечения
                    xf = fil.x;  yf = fil.y; zf = fil.z
                    cross_number_filler = index
                    break

        for index, sphere in enumerate(mol_list):
            cross_mol = 10 ** 4
            if index != cross_number_mol:
                D, a, b, c = is_collision(x0, y0, z0, sphere.x, sphere.y, sphere.z, xk, yk, zk, R)
                if D >= 0:
                    xk_m, yk_m, zk_m = collision(a, b, c, x0, y0, z0, xk, yk, zk)
                    cross_mol = math.sqrt((x0 - xk_m) ** 2 + (y0 - yk_m) ** 2 + (z0 - zk_m) ** 2)  #Расстояние от полости до точки пересечения
                    mol_volume = sphere.volume
                    xm = sphere.x; ym = sphere.y; zm = sphere.z
                    cross_number_mol = index
                    break


        if cross_filler < cross_mol:
            xk = xk_f; yk = yk_f; zk = zk_f
            file_arrow.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '   ')
            cos_alfa, cos_betta, cos_gamma = cos(x0, y0, z0, xk, yk, zk)
            H = S_for_new_koord(xk, yk, zk, x0, y0, z0, xf, yf, zf, R, r)
            old_x0, old_y0, old_z0 = old_koord(x0, y0, z0)
            x0, y0, z0 = new_koord(xk, yk, zk, H, cos_alfa, cos_betta, cos_gamma)
            Way += math.sqrt((x0 - old_x0) ** 2 + (y0 - old_y0) ** 2 + (z0 - old_z0) ** 2) #Расчет пути полости
            file_micropore.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
            file_arrow.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
            if capture == 1:
                file_capture_micropore.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
        elif cross_filler > cross_mol:
            xk = xk_m; yk = yk_m; zk = zk_m
            file_arrow.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '   ')
            if capture == 0:  #если молекула еще не захвачена
                if volume_cavity >= mol_volume:  #если объем полости больше объема молекулы
                    capture = 1  #Захват
                    old_x0, old_y0, old_z0 = old_koord(x0, y0, z0)
                    x0 = xm; y0 = ym; z0 = zm
                    Way += math.sqrt((x0 - old_x0) ** 2 + (y0 - old_y0) ** 2 + (z0 - old_z0) ** 2)  #Расчет пути полости
                    file_capture_micropore.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
                    mol_list.pop(cross_number_mol)
                    file_micropore.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
                    file_arrow.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
                else: #Не захват
                    cos_alfa, cos_betta, cos_gamma = cos(x0, y0, z0, xk, yk, zk)
                    H = S_for_new_koord(xk, yk, zk, x0, y0, z0, xf, yf, zf, mol_R, r)  #
                    old_x0, old_y0, old_z0 = old_koord(x0, y0, z0)
                    x0, y0, z0 = new_koord(xk, yk, zk, H, cos_alfa, cos_betta, cos_gamma)
                    Way += math.sqrt((x0 - old_x0) ** 2 + (y0 - old_y0) ** 2 + (z0 - old_z0) ** 2)  #Расчет пути полости
                    file_micropore.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
                    file_arrow.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
            else:  #если молекула уже захвачена
                cos_alfa, cos_betta, cos_gamma = cos(x0, y0, z0, xk, yk, zk)
                H = S_for_new_koord(xk, yk, zk, x0, y0, z0, xm, ym, zm, mol_R, r)
                old_x0, old_y0, old_z0 = old_koord(x0, y0, z0)
                x0, y0, z0 = new_koord(xk, yk, zk, H, cos_alfa, cos_betta, cos_gamma)
                Way += math.sqrt((x0 - old_x0) ** 2 + (y0 - old_y0) ** 2 + (z0 - old_z0) ** 2)  #Расчет пути полости
                file_micropore.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
                file_arrow.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
                file_capture_micropore.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '\n')
        elif zk > border and capture == 1:
            print ('Выход на поверхность')
            file_arrow.write(str(x0) + '   ' + str(y0) + '   ' + str(z0) + '   ' + str(xk) + '   ' + str(yk) + '   ' + str(zk))
            file_capture_micropore.write(str(xk) + '   ' + str(yk) + '   ' + str(zk))
            file_micropore.write(str(xk) + '   ' + str(yk) + '   ' + str(zk))
            break

    print(Way)
    file_micropore.close()
    file_filler.close()
    file_arrow.close()
    file_capture_micropore.close()
    file_molecules.close()


if __name__ == '__main__':
    main()




