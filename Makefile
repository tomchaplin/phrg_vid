.DEFAULT_GOAL := low

SCENE = IntroScene

low:
	manim -pql scene.py $(SCENE)

high:
	manim -pqh scene.py $(SCENE)

4k:
	manim -pqk scene.py $(SCENE)

play_low:
	mpv media/videos/scene/480p15/$(SCENE).mp4

play_high:
	mpv media/videos/scene/1080p60/$(SCENE).mp4

play_4k:
	mpv media/videos/scene/2160p60/$(SCENE).mp4

snap:
	manim -spql scene.py $(SCENE)
