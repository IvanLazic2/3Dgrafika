import math

file = open("cilindar.obj", "w")

numofdots = 18

r = 1
res = 360/numofdots

def genCircle(angle, z):
    if angle == 360:
        return

    x = r * math.cos(angle/360 * 2 * math.pi)
    y = r * math.sin(angle/360 * 2 * math.pi)

    file.write("v" + " " + str(x) + " " + str(y) + " " + str(z) + "\n")
    file.write("v" + " " + str(x) + " " + str(y) + " " + str(z) + "\n")

    genCircle(angle + res, z)


def genNorm(angle):
    if angle == 360:
        return

    file.write("vn " + str(math.cos(angle/360 * 2 * math.pi)) + " " + str(math.sin(angle/360 * 2 * math.pi)) + " " + "0\n")
    file.write("vn 0 0 1\n")

    genNorm(angle + res)


def genFaces(a, b, c, d, e, f, normal):
    if normal == numofdots * 2 + 1:
        return

    if normal == numofdots * 2 - 1:
        if a == numofdots * 2 - 1:
            a = 1
        if d == numofdots * 2 + 1:
            d = numofdots * 4 - 1
        if e == numofdots * 4 + 1:
            e = numofdots * 2 - 1
        if f == numofdots * 4 - 1:
            f = 1

    file.write("f " + str(a) + "//" + str(normal) + " " + str(b) + "//" + str(normal) + " " + str(c) + "//" + str(normal) + "\n")
    file.write("f " + str(d) + "//" + str(normal) + " " + str(e) + "//" + str(normal) + " " + str(f) + "//" + str(normal) + "\n")

    genFaces(a + 2, b + 2, c + 2, d + 2, e + 2, f + 2, normal + 2)


def genBase(a, b, normal, start):
    if normal == numofdots * 2:
        return

    file.write("f " + str(start) + "//" + str(normal) + " " + str(b) + "//" + str(normal) + " " + str(a) + "//" + str(normal) + "\n")

    genBase(a + 2, b + 2, normal + 2, start)


if __name__ == "__main__":
    genCircle(0, -1)
    genCircle(0, 1)

    genNorm(0)

    genFaces(1, 3, numofdots * 2 + 1, 3, numofdots * 2 + 3, numofdots * 2 + 1, 1)

    genBase(2, 4, 2, 2)
    genBase(numofdots * 2 + 4, numofdots * 2 + 2, 2, numofdots * 2 + 2)

    file.close()
