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


def getSquareGraph(directed=True, H0=2.5, W0=3.5):
    H = [2, 2, 0, 0]
    W = [0, 2, 2, 0]
    labels = ["a", "b", "c", "d"]
    V = [
        LabeledDot(labels[i]).shift(H0 * DOWN + W0 * RIGHT + H[i] * UP + W[i] * RIGHT)
        for i in range(4)
    ]
    A = [
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 1, 0, 1],
        [1, 0, 0, 0],
    ]
    EdgeFunc = Arrow if directed else Line
    E = [
        EdgeFunc(start=V[i], end=V[j], buff=0.1)
        for i in range(4)
        for j in range(4)
        if A[i][j] == 1
    ]
    return {"V": V, "E": E}


def getInitialGraph(directed=True, H0=1, W0=-1):
    H = [0, 2, 4, 4, 2, 0]
    W = [0, 0, 0, 2, 2, 2]
    V = [
        LabeledDot(f"v_{i}").shift(H0 * DOWN + W0 * RIGHT + H[i] * UP + W[i] * RIGHT)
        for i in range(6)
    ]
    A = [
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 1, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 1, 0, 0, 0, 0],
        [1, 0, 0, 0, 1, 0],
    ]
    EdgeFunc = Arrow if directed else Line
    E = [
        EdgeFunc(start=V[i], end=V[j], buff=0.1)
        for i in range(6)
        for j in range(6)
        if A[i][j] == 1
    ]
    smplxs = [
        (
            Polygon(*[v.get_center() for v in vs], color=BLUE)
            .set_stroke(width=0)
            .set_fill(color=BLUE, opacity=0.5)
        )
        for vs in [[V[1], V[2], V[3]], [V[1], V[4], V[3]]]
    ]
    return {"V": V, "E": E, "smplxs": smplxs}


class IntroScene(Scene):
    def construct(self):
        normal_graph = getInitialGraph(directed=False, H0=1, W0=-4)
        digraph = getInitialGraph(directed=True, H0=1, W0=-4)
        groups = [f"X_{i}" for i in range(CD_LEN)]
        digroups = [f"\overrightarrow{{X}}_{i}" for i in range(CD_LEN)]
        hom_groups = [f"H_{i}" for i in range(CD_LEN)]
        dihom_groups = [f"\overrightarrow{{H}}_{i}" for i in range(CD_LEN)]
        boundaries = [
            f"\\xrightarrow{{\\quad\\partial_{i+1}\\quad}}" for i in range(CD_LEN - 1)
        ]
        clique_comp = chainComplexFactory(groups, boundaries).shift(2.5 * DOWN)
        dflag_comp = chainComplexFactory(digroups, boundaries).shift(2.5 * DOWN)
        hom_comp = chainComplexFactory(hom_groups, boundaries).shift(2.5 * DOWN)
        dflag_hom_comp = chainComplexFactory(dihom_groups, boundaries).shift(2.5 * DOWN)
        clique_example = MathTex(r"\{v_1, v_2, v_3\}, \{v_1, v_3, v_4\}").shift(
            2 * RIGHT + 0.8 * DOWN
        )
        dclique_example = MathTex(r"(v_1, v_2, v_3)").shift(2 * RIGHT + 0.8 * DOWN)
        dclique_def1 = MathTex(r"(v_0, \dots , v_k)").shift(2 * RIGHT + 2 * UP)
        dclique_def2 = MathTex(r"i< j\ \implies\ v_i \to v_j").next_to(
            dclique_def1, DOWN
        )
        # Draw undirected graph
        self.play(*[Create(elem) for elem in normal_graph["V"] + normal_graph["E"]])
        self.wait(1)
        # Draw cliques as simplices
        self.bring_to_back(*normal_graph["smplxs"])
        self.play(*[Create(s) for s in normal_graph["smplxs"]], Write(clique_example))
        self.wait(1)
        # Chain complex
        self.play(*[Write(s) for s in clique_comp])
        self.wait(1)
        # Leads to homology
        self.play(
            *[Transform(clique_comp[2 * i], hom_comp[2 * i]) for i in range(CD_LEN)]
            + [
                Uncreate(clique_comp[2 * i + 1], run_time=0.5)
                for i in range(CD_LEN - 1)
            ]
        )
        self.wait(1)
        # Remove homology and simplices
        self.play(
            *[Uncreate(clique_comp[2 * i]) for i in range(CD_LEN)]
            + [Uncreate(s) for s in normal_graph["smplxs"]]
            + [Uncreate(clique_example)]
        )
        self.wait(1)
        # Make graph directed
        self.play(
            *[Transform(e1, e2) for e1, e2 in zip(normal_graph["E"], digraph["E"])]
        )
        self.wait(1)
        # Draw the directed clique
        self.play(Write(dclique_def1), Write(dclique_def2))
        self.wait(1)
        self.bring_to_back(digraph["smplxs"][0])
        self.play(Write(dclique_example), Create(digraph["smplxs"][0]))
        self.wait(1)
        # Draw the directed flag complex
        self.play(*[Write(s) for s in dflag_comp])
        self.wait(1)
        # Leads to homology
        self.play(
            *[
                Transform(dflag_comp[2 * i], dflag_hom_comp[2 * i])
                for i in range(CD_LEN)
            ]
            + [Uncreate(dflag_comp[2 * i + 1], run_time=0.5) for i in range(CD_LEN - 1)]
        )
        self.wait(1)


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
            Tex(r"$abdc$\\", r"is ", r"not", r" allowed")
            .set_color_by_tex("not", RED)
            .next_to(allowed_path, DOWN)
            .align_to(allowed_path, LEFT)
            .shift(1.5 * DOWN + 2 * RIGHT)
        )
        square_graph = getSquareGraph()
        sg_green_edges = [0, 1]
        sg_red_edges = [3]
        all_defs = MathTex(
            r"\Lambda_p &:= \mathbb{Z}\langle\{\text{elementary }p\text{-paths}\}\rangle\\",
            r"\mathcal{A}_p &:= \mathbb{Z}\langle\{\text{allowed }p\text{-paths}\}\rangle\\",
            r"\partial_p(v_0 \dots v_p) &:= \sum_{i=0}^p (-1)^i v_0 \dots \hat{v_i} \dots v_p",
        ).shift(2 * DOWN)

        self.add(graph_def, elem_path, allowed_path)
        self.wait(5)
        self.play(*[Create(elem) for elem in square_graph["V"] + square_graph["E"]])
        self.wait(1)
        self.play(
            *(
                [Write(not_allowed_eg)]
                + [
                    square_graph["E"][i].animate.set_fill(GREEN).set_stroke(GREEN)
                    for i in sg_green_edges
                ]
                + [
                    square_graph["E"][i].animate.set_fill(RED).set_stroke(RED)
                    for i in sg_red_edges
                ]
            )
        )
        self.wait(3)
        self.play(
            *[
                Uncreate(elem)
                for elem in [not_allowed_eg] + square_graph["V"] + square_graph["E"]
            ]
        )
        self.play(Write(all_defs))
        self.wait(5)


class OmegaDef(Scene):
    def construct(self):
        lambdas = [f"\\Lambda_{i}" for i in range(CD_LEN)]
        As = [f"\\mathcal{{A}}_{i}" for i in range(CD_LEN)]
        Omegas = [f"\\Omega_{i}" for i in range(CD_LEN)]
        Homs = [f"H^{{\Xi}}_{i}" for i in range(CD_LEN)]
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
        Omegas_complex_centred = chainComplexFactory(Omegas, boundaries).align_to(
            lambda_complex, LEFT
        )
        Hom_complex = chainComplexFactory(Homs, boundaries).align_to(
            lambda_complex, LEFT
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
        self.play(
            *[
                Transform(As_complex[2 * i + 1], Omegas_complex_G[2 * i + 1])
                for i in range(CD_LEN - 1)
            ],
            Transform(a2box, o1goodbox),
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
        # Centre Omega and delete everything else
        self.play(
            *[FadeOut(i) for i in ls_group + hooks],
            FadeOut(ls_problem),
            FadeOut(lambda_complex),
            FadeOut(a2box),
        )
        self.play(
            *[
                Transform(As_complex[2 * i + 1], Omegas_complex_centred[2 * i + 1])
                for i in range(CD_LEN - 1)
            ]
            + [
                ApplyMethod(As_complex[2 * i].match_y, Omegas_complex_centred[2 * i])
                for i in range(CD_LEN)
            ]
        )
        self.wait(1)
        self.play(
            *[Uncreate(As_complex[2 * i + 1]) for i in range(CD_LEN - 1)]
            + [Transform(As_complex[2 * i], Hom_complex[2 * i]) for i in range(CD_LEN)]
        )
        self.wait(1)
