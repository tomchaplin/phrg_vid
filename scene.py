from manim import *

CD_LEN = 4


def my_interlace(l1, l2):
    return list(sum(zip(l1, l2 + [0]), ())[:-1])


def chainComplexFactory(groups, maps, group_width="25pt", complex_type="chain"):
    all_text = my_interlace(
        [f"\\parbox{{ {group_width} }}{{\\centering ${g}$ }}" for g in groups],
        [f"${m}$" for m in maps],
    )
    if complex_type != "cochain":
        all_text.reverse()
    return Tex(*all_text)


def getHookAbove(thingBelow):
    inclds = MathTex(r"\lhook\joinrel\longrightarrow")
    inclds.rotate(3 * PI / 2)
    inclds.next_to(thingBelow, UP)
    return inclds


def getRedArrowLeftOf(thingRight, label):
    arr = MathTex(label, color=RED)
    arr.next_to(thingRight, LEFT, buff=0.1)
    return arr


def getLongSquare():
    a = Dot(DOWN * 3 + LEFT * 5)
    b = Dot(DOWN * 2 + LEFT * 5)
    c = Dot(DOWN * 3 + LEFT * 4)
    d = Dot(DOWN * 2 + LEFT * 4)
    a_label = MathTex("a").next_to(a, LEFT).scale(0.7)
    b_label = MathTex("b").next_to(b, LEFT).scale(0.7)
    c_label = MathTex("c").next_to(c, RIGHT).scale(0.7)
    d_label = MathTex("d").next_to(d, RIGHT).scale(0.7)
    Aab = Arrow(start=a, end=b, buff=0.1)
    Aac = Arrow(start=a, end=c, buff=0.1)
    Abd = Arrow(start=b, end=d, buff=0.1)
    Acd = Arrow(start=c, end=d, buff=0.1)
    return [a, b, c, d, Aab, Aac, Abd, Acd, a_label, b_label, c_label, d_label]


def getLongSquareProblem():
    return (
        MathTex(r"\partial(abd) = bd - ", r"ad", r" + ab")
        .set_color_by_tex("ad", RED)
        .shift(DOWN * 2.5 + RIGHT)
    )


def getLongSquareSolution():
    return (
        MathTex(
            r"\partial(abd - acd) &= (bd - ",
            r"{ad}",
            r" + ab) \\&\quad\  - (cd - ",
            r"{ad}",
            r" + ac) \\ &= ab + bd - cd - ac",
        )
        .set_color_by_tex("{ad}", RED)
        .shift(DOWN * 2.5 + RIGHT)
    )


def getInterestingGraph():
    H = [2, 2, 0, 0]
    W = [0, 2, 2, 0]
    V = [
        LabeledDot(f"v_{i}").shift(2.5 * DOWN + 3.5 * RIGHT + H[i] * UP + W[i] * RIGHT)
        for i in range(4)
    ]
    A = [
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 1, 0, 1],
        [1, 0, 0, 0],
    ]
    E = [
        Arrow(start=V[i], end=V[j], buff=0.1)
        for i in range(4)
        for j in range(4)
        if A[i][j] == 1
    ]
    return {"V": V, "E": E}


class GroupDef(Scene):
    def construct(self):
        graph_def = Tex(
            r"Directed graph $G=(V, E)$ i.e.\ $E\subseteq V\times V$"
        ).to_corner(UL)
        elem_path = (
            Tex(
                "An ",
                r"elementary $p$-path ",
                r"is an ordered tuple of $(p+1)$ vertices\\$v_0\dots v_p$",
            )
            .set_color_by_tex("elementary", PURPLE)
            .next_to(graph_def, DOWN)
            .align_to(graph_def, LEFT)
        )
        allowed_path = (
            Tex(
                "An ",
                r"allowed $p$-path ",
                r"must follow the directed edges\\",
                r"i.e.\ $(v_{i-1},v_i)\in E$",
            )
            .set_color_by_tex("allowed", PURPLE)
            .next_to(elem_path, DOWN)
            .align_to(elem_path, LEFT)
        )
        not_allowed_eg = (
            Tex(r"$v_0 v_1 v_3 v_2$\\", r"is ", r"not", r" allowed")
            .set_color_by_tex("not", RED)
            .next_to(allowed_path, DOWN)
            .align_to(allowed_path, LEFT)
            .shift(1.5 * DOWN + 2 * RIGHT)
        )
        interesting_graph = getInterestingGraph()
        ig_green_edges = [0, 1]
        ig_red_edges = [3]
        all_defs = MathTex(
            r"\Lambda_p &:= R\langle\{\text{elementary }p\text{-paths}\}\rangle\\",
            r"\mathcal{A}_p &:= R\langle\{\text{allowed }p\text{-paths}\}\rangle\\",
            r"\partial_p(v_0 \dots v_p) &:= \sum_{i=0}^p (-1)^i v_0 \dots \hat{v_i} \dots v_p",
        ).shift(2 * DOWN)

        self.add(graph_def, elem_path, allowed_path)
        self.wait(5)
        self.play(
            *[Create(elem) for elem in interesting_graph["V"] + interesting_graph["E"]]
        )
        self.wait(1)
        self.play(
            *(
                [Write(not_allowed_eg)]
                + [
                    interesting_graph["E"][i].animate.set_fill(GREEN).set_stroke(GREEN)
                    for i in ig_green_edges
                ]
                + [
                    interesting_graph["E"][i].animate.set_fill(RED).set_stroke(RED)
                    for i in ig_red_edges
                ]
            )
        )
        self.wait(3)
        self.play(
            *[
                Uncreate(elem)
                for elem in [not_allowed_eg]
                + interesting_graph["V"]
                + interesting_graph["E"]
            ]
        )
        self.play(Write(all_defs))
        self.wait(5)


class OmegaDef(Scene):
    def construct(self):
        lambdas = [f"\\Lambda_{i}" for i in range(CD_LEN)]
        As = [f"\\mathcal{{A}}_{i}" for i in range(CD_LEN)]
        Omegas = [f"\\Omega_{i}" for i in range(CD_LEN)]
        boundaries = [
            f"\\xrightarrow{{\\quad\\partial_{i+1}\\quad}}" for i in range(CD_LEN - 1)
        ]

        lambda_complex = chainComplexFactory(lambdas, boundaries)
        hooks = [getHookAbove(lambda_complex[2 * i]) for i in range(CD_LEN)]
        As_complex = (
            chainComplexFactory(As, boundaries)
            .set_color_by_tex("xrightarrow", RED)
            .next_to(hooks[0], UP)
            .align_to(lambda_complex, LEFT)
        )
        Omegas_complex_R = (
            chainComplexFactory(Omegas, boundaries)
            .set_color_by_tex("xrightarrow", RED)
            .next_to(hooks[0], UP)
            .align_to(lambda_complex, LEFT)
        )
        Omegas_complex_G = (
            chainComplexFactory(Omegas, boundaries)
            .set_color_by_tex("xrightarrow", GREEN)
            .next_to(hooks[0], UP)
            .align_to(lambda_complex, LEFT)
        )

        a2box = SurroundingRectangle(As_complex[2], buff=0.1)
        l2box = SurroundingRectangle(lambda_complex[2], buff=0.1)
        l1box = SurroundingRectangle(lambda_complex[4], buff=0.1)
        a1boxbad = SurroundingRectangle(lambda_complex[4], color=RED, buff=0.1).shift(
            UP * 0.1
        )
        o1goodbox = SurroundingRectangle(Omegas_complex_G[4], color=GREEN, buff=0.1)

        Omega_def = MathTex(
            r"\Omega_p := \{v \in \mathcal{A}_p \ |\ \partial_p v \in \mathcal{A}_{p-1}\}"
        ).shift(DOWN * 2.5)

        # Draw initial Lambda Complex
        self.play(Write(lambda_complex))
        self.wait()
        # Show allowed complex included
        self.play(
            *[Write(h) for h in hooks]
            + [Write(As_complex[2 * i]) for i in range(CD_LEN)]
        )
        self.wait(2)
        # Show problem with boundary in allowed complex
        self.play(Create(a2box))
        self.wait()
        self.play(Transform(a2box, l2box))
        self.wait()
        self.play(Transform(a2box, l1box))
        self.wait()
        for i in range(2):
            self.play(Transform(a2box, a1boxbad), run_time=0.5)
            self.play(Transform(a2box, l1box), run_time=0.5)
        self.play(*[Write(As_complex[2 * i + 1]) for i in range(CD_LEN - 1)])
        self.wait(2)
        # Example 1
        ls_group = getLongSquare()
        ls_problem = getLongSquareProblem()
        self.play(*[Create(i) for i in ls_group])
        self.wait(2)
        self.play(Write(ls_problem))
        self.wait(5)
        self.play(*[Uncreate(i) for i in ls_group + [ls_problem]])
        # Define Omega complex
        self.play(Write(Omega_def))
        self.wait(5)
        self.play(Transform(As_complex, Omegas_complex_R))
        self.wait(2)
        self.play(Transform(a2box, o1goodbox))
        self.wait()
        self.play(
            *[
                Transform(As_complex[2 * i + 1], Omegas_complex_G[2 * i + 1])
                for i in range(CD_LEN - 1)
            ]
        )
        self.wait(5)
        self.play(Uncreate(Omega_def))
        # Example 1 fixed
        ls_group = getLongSquare()
        ls_problem = getLongSquareProblem()
        ls_sol = getLongSquareSolution()
        self.play(*[Create(i) for i in ls_group] + [Write(ls_problem)])
        self.wait(2)
        self.play(Transform(ls_problem, ls_sol))
        self.wait(5)
