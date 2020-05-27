from django.http import HttpResponse


def index(request):
    return HttpResponse("Starting ska sensitivity database application")
    