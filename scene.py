from manim import *
import copy

CD_LEN = 4


def fadeOutAll(self, run_time=0.5):
    return [FadeOut(mob, run_time=run_time) for mob in self.mobjects]


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


def getGrigoryanCite(paper="homology"):
    if paper == "homology":
        str = "(Grigor'yan et al. 2013) Homologies of path complexes and digraph"
    elif paper == "homotopy":
        str = "(Grigor'yan et al. 2014) Homotopy theory for digraphs"

    return Text(str).scale(0.3).to_corner(DL)


def getLongSquare(offset=0 * DOWN):
    a = Dot(offset + DOWN * 2.5 + LEFT * 5)
    b = Dot(offset + DOWN * 1.5 + LEFT * 5)
    c = Dot(offset + DOWN * 2.5 + LEFT * 4)
    d = Dot(offset + DOWN * 1.5 + LEFT * 4)
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
        .shift(DOWN * 2 + RIGHT)
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
        .shift(DOWN * 2 + RIGHT)
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


def getBigCycle():
    V = [
        Dot(2 * (LEFT + UP)),
        Dot(2 * (RIGHT + UP)),
        Dot(3 * RIGHT),
        Dot(2 * (RIGHT + DOWN)),
        Dot(2 * (LEFT + DOWN)),
        Dot(3 * LEFT),
    ]
    flip = [True, True, False, True, False]
    E = [
        Arrow(V[i], V[i + 1]) if flip[i] else Arrow(V[i + 1], V[i]) for i in range(5)
    ] + [Arrow(V[5], V[0])]
    dirs = [UP, UP, RIGHT, DOWN, DOWN, LEFT]
    labels = [MathTex(f"v_{i}").next_to(V[i], dirs[i]) for i in range(len(V))]
    return {"V": V, "E": E, "labels": labels}


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
        # Fade all out
        self.play(
            *[FadeOut(e, run_time=0.5) for e in dflag_comp]
            + [
                Uncreate(e, run_time=0.5)
                for e in (
                    normal_graph["V"] + normal_graph["E"] + normal_graph["smplxs"]
                )
            ]
            + [
                FadeOut(e, run_time=0.5)
                for e in [dclique_example, dclique_def1, dclique_def2]
            ]
        )


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
            .scale(0.8)
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
            .scale(0.8)
            .next_to(allowed_path, DOWN)
            .align_to(allowed_path, LEFT)
            .shift(1.5 * DOWN + 2 * RIGHT)
        )
        citation = getGrigoryanCite()
        square_graph = getSquareGraph()
        sg_green_edges = [0, 1]
        sg_red_edges = [3]
        all_defs = (
            MathTex(
                r"\Lambda_p &:= \mathbb{Z}\langle\{\text{elementary }p\text{-paths}\}\rangle\\",
                r"\mathcal{A}_p &:= \mathbb{Z}\langle\{\text{allowed }p\text{-paths}\}\rangle\\",
                r"\partial_p(v_0 \dots v_p) &:= \sum_{i=0}^p (-1)^i v_0 \dots \hat{v_i} \dots v_p",
            )
            .shift(1.1 * DOWN)
            .scale(0.8)
        )

        self.play(
            *[
                FadeIn(e, run_time=0.5)
                for e in [graph_def, elem_path, allowed_path, citation]
            ]
        )
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
                Uncreate(elem, run_time=0.5)
                for elem in [not_allowed_eg] + square_graph["V"] + square_graph["E"]
            ]
        )
        self.play(FadeIn(all_defs))
        self.wait(5)
        self.play(
            *[
                FadeOut(e, run_time=0.5)
                for e in [all_defs, graph_def, elem_path, allowed_path]
            ]
        )


class OmegaDef(Scene):
    def construct(self):
        lambdas = [f"\\Lambda_{i}" for i in range(CD_LEN)]
        As = [f"\\mathcal{{A}}_{i}" for i in range(CD_LEN)]
        Omegas = [f"\\Omega_{i}" for i in range(CD_LEN)]
        Homs = [f"H^{{\Xi}}_{i}" for i in range(CD_LEN)]
        boundaries = [
            f"\\xrightarrow{{\\quad\\partial_{i+1}\\quad}}" for i in range(CD_LEN - 1)
        ]

        citation = getGrigoryanCite()
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
        ).shift(DOWN * 2)

        # Draw initial Lambda Complex
        self.add(citation)
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
        self.play(FadeOut(Omega_def, run_time=0.5))
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
        # Fade out hom groups
        self.play(*fadeOutAll(self))


class Generators(Scene):
    def construct(self):
        citation = getGrigoryanCite(paper="homotopy")
        first_two = (
            MathTex(
                r"\Omega_0 &= \mathbb{Z}\langle V \rangle\\",
                r"\Omega_1 &= \mathbb{Z}\langle E \rangle",
            )
            .to_corner(UL)
            .shift(DOWN + 0.5 * RIGHT)
        )
        omega_2 = MathTex("\Omega_2 = \mathbb{Z}\Bigg\langle").align_to(first_two, LEFT)
        eyeglass = {}
        eyeglass["V"] = [
            Dot(3.7 * LEFT),
            Dot(2.7 * LEFT),
            Dot(1.7 * LEFT),
        ]
        eyeglass["E"] = [
            e.scale(0.6)
            for e in [
                CurvedArrow(
                    eyeglass["V"][0].get_center(),
                    eyeglass["V"][1].get_center(),
                ),
                CurvedArrow(
                    eyeglass["V"][1].get_center(),
                    eyeglass["V"][0].get_center(),
                ),
                CurvedArrow(
                    eyeglass["V"][1].get_center(),
                    eyeglass["V"][2].get_center(),
                ),
                CurvedArrow(
                    eyeglass["V"][2].get_center(),
                    eyeglass["V"][1].get_center(),
                ),
            ]
        ]
        names = ["i", "j", "k"]
        eyeglass["labels"] = [
            MathTex(names[i]).next_to(eyeglass["V"][i], UP).scale(0.7) for i in range(3)
        ]
        eyeglass["element"] = (
            MathTex("jij - jkj").next_to(eyeglass["V"][1], 1.2 * DOWN).scale(0.7)
        )

        triangle = {}
        triangle["V"] = [Dot(0.2 * LEFT), Dot(0.55 * RIGHT + UP), Dot(1.3 * RIGHT)]
        triangle["E"] = [
            Arrow(start=triangle["V"][0], end=triangle["V"][1]),
            Arrow(start=triangle["V"][0], end=triangle["V"][2]),
            Arrow(start=triangle["V"][1], end=triangle["V"][2]),
        ]
        triangle["labels"] = [
            MathTex(names[i]).next_to(triangle["V"][i], [LEFT, UP, RIGHT][i]).scale(0.7)
            for i in range(3)
        ]
        triangle["element"] = MathTex("ijk").next_to(triangle["E"][1], DOWN).scale(0.7)
        ls = getLongSquare(offset=8.7 * RIGHT + 2.5 * UP)
        ls_label = MathTex("abd - acd").next_to(ls[5], DOWN).scale(0.7)

        omega_2close = MathTex(r"\Bigg\rangle").next_to(omega_2, 40 * RIGHT)

        self.play(Write(first_two), FadeIn(citation))
        self.wait(1)
        self.play(Write(omega_2))
        for gen in [eyeglass, triangle]:
            self.play(
                *[
                    Create(e)
                    for e in gen["V"] + gen["E"] + gen["labels"] + [gen["element"]]
                ]
            )
            self.wait(1)

        self.play(*[Create(e) for e in ls + [ls_label, omega_2close]])
        self.wait(5)
        self.play(*fadeOutAll(self))


class ERGraph(Scene):
    def construct(self):
        V = [
            Dot(DOWN + 2 * (LEFT + UP)),
            Dot(DOWN + 2 * (RIGHT + UP)),
            Dot(DOWN + 2 * (RIGHT + DOWN)),
            Dot(DOWN + 2 * (LEFT + DOWN)),
            Dot(DOWN + 3 * LEFT),
            Dot(DOWN + 3 * RIGHT),
        ]
        arrows = [
            Arrow(2 * (LEFT + DOWN) + DOWN, 3 * RIGHT + DOWN, color=MAROON_A),
            Arrow(2 * (LEFT + UP) + DOWN, 2 * (LEFT + DOWN) + DOWN, color=MAROON_A),
            Arrow(3 * LEFT + DOWN, 2 * (LEFT + UP) + DOWN, color=MAROON_A),
            CurvedArrow(2 * (LEFT + UP) + DOWN, 3 * LEFT + DOWN, color=MAROON_A),
        ]
        name = MathTex(r"\overrightarrow{G}(n,p)").to_corner(UL)
        nodes = Tex("$n$ nodes").next_to(name, DOWN).align_to(name, LEFT)
        prob = (
            Tex("Edges appear indepdently with probability $p$", color=MAROON_A)
            .next_to(nodes, DOWN)
            .align_to(nodes, LEFT)
        )
        self.play(FadeIn(name))
        self.play(*[FadeIn(v) for v in V], Write(nodes))
        for i in range(len(arrows)):
            self.play(Create(arrows[i]), run_time=0.5)
        self.play(Write(prob, run_time=1))
        self.wait(1)
        self.play(*fadeOutAll(self))


class EmpiricalDist(Scene):
    def construct(self):
        im = ImageMobject("assets/prob_nonzero.png")
        im.scale(0.3)
        l1 = Line(
            start=(1.5 * UP + 3.1 * LEFT),
            end=(2 * RIGHT),
            color=RED,
            stroke_width=15,
        ).set_opacity(0.9)
        l2 = Line(
            start=(3 * UP + 3.1 * LEFT),
            end=(2 * RIGHT + 2.2 * UP),
            color=RED,
            stroke_width=15,
        ).set_opacity(0.9)
        l3 = Line(
            start=(2.8 * UP + 3.1 * LEFT), end=(2.8 * UP + 2 * RIGHT), color=BLACK
        )
        l4 = Line(
            start=(2.8 * UP + 3.1 * LEFT), end=(1.5 * UP + 2 * RIGHT), color=BLACK
        )
        l5 = Line(
            start=(2.8 * UP + 3.1 * LEFT), end=(3 * DOWN + 2 * RIGHT), color=BLACK
        )
        pt1 = Dot(DOWN + LEFT, radius=0.1, color=RED)
        pt2 = Dot(2.9 * UP + RIGHT, radius=0.1, color=RED)
        pt3 = Dot(1.6 * UP + 0.6 * RIGHT, radius=0.1, color=RED)
        self.play(FadeIn(im), run_time=0.5)
        self.play(Create(l1), Create(l2))
        self.wait(2)
        self.play(Create(pt1))
        self.wait(3)
        self.play(Transform(pt1, pt2))
        self.wait(3)
        self.play(Transform(pt1, pt3))
        self.wait(3)
        self.play(Uncreate(pt1))
        self.play(Create(l3))
        self.wait(1)
        self.play(Transform(l3, l4))
        self.wait(1)
        self.play(Transform(l3, l5))
        self.wait(3)
        self.play(*fadeOutAll(self))


class SmallDensities(Scene):
    def construct(self):
        bigCycle = getBigCycle()
        prob_nocycle = (
            MathTex(
                r"&\mathbb{P}(\overrightarrow{\beta}_1(G)>0)\\",
                r"&\quad\leq\mathbb{P}(\exists\text{ undirected cycle})\\",
                r"&\quad\leq\sum_{L=2}^\infty \binom{n}{L}L! 2^L p^L",
                r"\leq\frac{(2np)^2}{1-2np}",
                r"\to 0",
            )
            .to_corner(UL)
            .shift(DOWN + RIGHT)
        )
        p_condition = (
            Tex(r"So long as $p=n^\alpha$, with $\alpha > -1$")
            .next_to(prob_nocycle, DOWN)
            .align_to(prob_nocycle, LEFT)
        )
        gradient_implication = (
            Tex(r"$\therefore$ The gradient of the top boundary is $\geq -1$")
            .next_to(p_condition, DOWN)
            .align_to(p_condition, LEFT)
        )
        self.play(*[Create(v) for v in bigCycle["V"]])
        for i in range(len(bigCycle["E"])):
            self.play(Create(bigCycle["E"][i], run_time=0.5))
        self.wait(1)
        self.play(*fadeOutAll(self, run_time=1))
        for i in range(len(prob_nocycle)):
            self.play(Write(prob_nocycle[i]))
            self.wait(1)
        self.play(Write(p_condition))
        self.wait(1)
        self.play(*fadeOutAll(self))
        # self.wait(1)
        # self.play(Write(gradient_implication))
        # self.wait(1)


class LargeDensities(Scene):
    def construct(self):
        bigCycle = getBigCycle()
        centre = Dot(ORIGIN)
        centrep = Dot(5.255 * RIGHT + 3.25 * UP)
        centre_edges = [
            Arrow(bigCycle["V"][0], centre),
            Arrow(centre, bigCycle["V"][2]),
            Arrow(centre, bigCycle["V"][3]),
        ]
        centrep_edges = [
            Arrow(bigCycle["V"][0], centrep).set_color(GREEN),
            Arrow(centrep, bigCycle["V"][2]),
            Arrow(centrep, bigCycle["V"][3]).set_color(GREEN),
        ]
        homologous_stmt = (
            MathTex("v_0 v_1 + v_1 v_2 - v_3 v_2", r"\sim", r"v_0\kappa + \kappa v_3")
            .set_color_by_tex("v_1", RED)
            .set_color_by_tex("kappa", GREEN)
        ).to_edge(DOWN)
        homologous_stmtp = (
            MathTex("v_0 v_1 + v_1 v_2 - v_3 v_2", r"\sim", r"v_0\kappa' + \kappa' v_3")
            .set_color_by_tex("v_1", RED)
            .set_color_by_tex("kappa", GREEN)
        ).to_edge(DOWN)
        kappa = MathTex(r"\kappa").next_to(centre, DOWN + LEFT)
        kappap = MathTex(r"\kappa'").next_to(centrep, RIGHT)

        prob_nodcentre = MathTex(
            r"&\mathbb{P}(\exists\text{ undirected 3-path without a directed centre})\\",
            r"&\quad\leq 4n^4 p^3 e^{-p^3(n-4)}",
            r"\to 0",
        )
        p_condition = (
            Tex(r"So long as $p=n^\alpha$, with $\alpha < -1/3$")
            .next_to(prob_nodcentre, DOWN)
            .align_to(prob_nodcentre, LEFT)
        )
        gradient_implication = (
            Tex(r"$\therefore$ The gradient of the top boundary is $\leq -1/3$")
            .next_to(p_condition, DOWN)
            .align_to(p_condition, LEFT)
        )

        self.play(
            *[FadeIn(e) for e in bigCycle["V"] + bigCycle["E"] + bigCycle["labels"]],
        )
        self.wait(1)
        self.play(Create(centre), Write(kappa))
        self.wait(1)
        self.play(*[Create(e) for e in centre_edges])
        self.wait(1)
        # Highlight subgraph by lowering opacity of other edges
        self.play(*[e.animate.set_opacity(0.2) for e in bigCycle["E"][3:]])
        self.wait(0.5)
        self.play(*[e.animate.set_opacity(1) for e in bigCycle["E"][3:]])
        self.wait(1)
        self.play(
            *[
                Transform(e, e.set_color(RED))
                for e in [bigCycle["E"][i] for i in range(3)]
            ]
            + [
                Transform(centre_edges[i], centre_edges[i].set_color(GREEN))
                for i in [0, 2]
            ],
            Write(homologous_stmt),
        )
        self.wait(3)
        self.play(
            *[Transform(e1, e2) for e1, e2 in zip(centre_edges, centrep_edges)],
            Transform(centre, centrep),
            Transform(kappa, kappap),
            Transform(homologous_stmt[2], homologous_stmtp[2]),
        )
        self.wait(1)
        self.play(*fadeOutAll(self, run_time=1))
        self.wait(1)
        for i in range(len(prob_nodcentre)):
            self.play(Write(prob_nodcentre[i]))
            if i != 0:
                self.wait(1)
        self.play(Write(p_condition))
        self.wait(1)
        self.play(*fadeOutAll(self))
        # self.wait(1)
        # self.play(Write(gradient_implication))
        # self.wait(1)


class CycleReduction(Scene):
    def construct(self):
        bigCycle = getBigCycle()
        extra_V = [Dot(RIGHT), Dot(DOWN), Dot(1.5 * LEFT + 0.5 * DOWN)]
        extra_E = [
            Arrow(bigCycle["V"][0], extra_V[0]),
            Arrow(extra_V[0], bigCycle["V"][3]),
            Arrow(bigCycle["V"][0], extra_V[1]),
            Arrow(extra_V[1], bigCycle["V"][4]),
            Arrow(bigCycle["V"][0], extra_V[2]),
            Arrow(extra_V[2], bigCycle["V"][5]),
        ]
        centred_V = [Dot(2 * UP), Dot(DOWN + 2 * RIGHT), Dot(DOWN + 2 * LEFT)]
        centred_E = [
            Arrow(centred_V[0], centred_V[1]),
            Arrow(centred_V[1], centred_V[2]),
            Arrow(centred_V[2], centred_V[0]),
        ]
        final_V = [bigCycle["V"][0], extra_V[2], bigCycle["V"][5]]
        final_E = [extra_E[4], extra_E[5], bigCycle["E"][5]]
        nul_E = Arrow(ORIGIN, ORIGIN).scale(0)
        hom0_V = Dot(ORIGIN)
        hom0_E = [
            DashedLine(hom0_V, centred_V[i]).scale(0.8).add_tip() for i in range(3)
        ]
        self.play(
            *[FadeIn(e) for e in bigCycle["V"] + bigCycle["E"]],
        )
        self.wait(1)
        self.play(*[FadeIn(extra_V[0])] + [FadeIn(extra_E[i]) for i in range(2)])
        self.play(
            *[FadeOut(bigCycle["E"][i]) for i in range(3)]
            + [FadeOut(bigCycle["V"][i + 1]) for i in range(2)]
        )
        self.play(*[FadeIn(extra_E[i + 2]) for i in range(2)], FadeIn(extra_V[1]))
        self.play(
            FadeOut(bigCycle["V"][3]),
            FadeOut(extra_V[0]),
            FadeOut(bigCycle["E"][3]),
            *[FadeOut(extra_E[i]) for i in range(2)],
        )
        self.play(*[FadeIn(extra_E[i + 4]) for i in range(2)], FadeIn(extra_V[2]))
        self.play(
            FadeOut(bigCycle["V"][4]),
            FadeOut(extra_V[1]),
            FadeOut(bigCycle["E"][4]),
            *[FadeOut(extra_E[i + 2]) for i in range(2)],
        )
        self.play(
            *[Transform(e1, e2) for e1, e2 in zip(final_E, centred_E)]
            + [Transform(v1, v2) for v1, v2 in zip(final_V, centred_V)]
        )
        self.play(FadeIn(hom0_V), *[FadeIn(e) for e in hom0_E])
        self.play(
            *[Transform(e, nul_E) for e in hom0_E + final_E]
            + [Transform(v, hom0_V) for v in final_V]
        )
        self.wait(1)
        self.play(*fadeOutAll(self))
