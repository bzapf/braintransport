from ufl import exp, ln


# def sigmoid(x):
#     try:
#         retval = 1 / (1 + exp(-x))
#     except:
#         x = x[0]
#         retval = [1 / (1 + exp(-x))]
#     return retval


# def inverse_sigmoid(y):
#     try:
#         retval = -1 * ln((1 / y) - 1)
#     except:
#         y = y[0]
#         retval = [-1 * ln((1 / y) - 1)]
#     return retval    


class PositiveParam():
    def __init__(self, min=0.1):
        self.min = min



    def __call__(self, delta):           
        return self.f(delta)


    def f(self, delta):
        return delta

    def f_inverse(self, param):
        return param

    # def f(self, delta):
    #     return (delta) ** 2 + self.min

    # def f_inverse(self, param):
    #     return (param - self.min) ** (1 / 2)



if __name__ == "__main__":
    d0 = 1
    d = PositiveParam(d0)
    print("D=", d0)
    print("delta = f_inverse(d0)", PositiveParam.f_inverse(d0))
    print("D(delta)=", PositiveParam.f(PositiveParam.f_inverse(d0)), "(should be ", d0, ")")