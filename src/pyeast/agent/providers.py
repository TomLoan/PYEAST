"""LLM provider abstraction for the PYEAST agent.

Both providers exchange messages in Anthropic's content-block shape:
    {"role": "user" | "assistant", "content": [...]}
with blocks of type "text", "tool_use", or "tool_result". This is the
canonical format used by agent.py's conversation loop; OllamaProvider
translates to/from it at its own boundary.
"""

import json
import os
import sys
import uuid

import urllib.error
import urllib.request


class LLMProvider:
    """Interface for an LLM backend used by the PYEAST agent."""

    def send(self, messages: list, system_prompt: str, tools: list) -> dict:
        """Send the conversation so far and return a normalized response.

        Returns {"content": [...], "stop_reason": "end_turn" | "tool_use"}
        using Anthropic-style content blocks. Implementations are
        responsible for printing any text output themselves.
        """
        raise NotImplementedError


class AnthropicProvider(LLMProvider):
    def __init__(self, model: str):
        import anthropic

        api_key = os.environ.get("ANTHROPIC_API_KEY")
        if not api_key:
            print(
                "Error: ANTHROPIC_API_KEY not set.\n"
                "Add it to a .env file or export it as an environment variable.\n"
                "See .env.example for the expected format."
            )
            sys.exit(1)

        self.client = anthropic.Anthropic(api_key=api_key)
        self.model = model

    def send(self, messages: list, system_prompt: str, tools: list) -> dict:
        with self.client.messages.stream(
            model=self.model,
            max_tokens=32000,
            system=[{
                "type": "text",
                "text": system_prompt,
                "cache_control": {"type": "ephemeral"},
            }],
            tools=tools,
            messages=messages,
        ) as stream:
            for text in stream.text_stream:
                print(text, end="", flush=True)
            response = stream.get_final_message()

        content_blocks = []
        for block in response.content:
            if block.type == "text":
                content_blocks.append({"type": "text", "text": block.text})
            elif block.type == "tool_use":
                content_blocks.append({
                    "type": "tool_use",
                    "id": block.id,
                    "name": block.name,
                    "input": block.input,
                })

        return {"content": content_blocks, "stop_reason": response.stop_reason}


_KNOWN_TOOL_MODELS = "llama3.1, qwen2.5, mistral-nemo"


class OllamaProvider(LLMProvider):
    def __init__(self, model: str | None):
        if not model:
            print(
                "Error: --model is required for --provider ollama.\n"
                f"Use a tool-calling-capable model, e.g.: {_KNOWN_TOOL_MODELS}\n"
                "Pull one first with: ollama pull <model>"
            )
            sys.exit(1)

        self.model = model
        self.host = os.environ.get("OLLAMA_HOST", "http://localhost:11434").rstrip("/")

    @staticmethod
    def _to_ollama_tools(tools: list) -> list:
        return [
            {
                "type": "function",
                "function": {
                    "name": tool["name"],
                    "description": tool["description"],
                    "parameters": tool["input_schema"],
                },
            }
            for tool in tools
        ]

    @staticmethod
    def _to_ollama_messages(messages: list, system_prompt: str) -> list:
        ollama_messages = [{"role": "system", "content": system_prompt}]

        for message in messages:
            role = message["role"]
            content = message["content"]

            if isinstance(content, str):
                ollama_messages.append({"role": role, "content": content})
                continue

            text_parts = []
            tool_calls = []
            for block in content:
                if block["type"] == "text":
                    text_parts.append(block["text"])
                elif block["type"] == "tool_use":
                    tool_calls.append({
                        "function": {"name": block["name"], "arguments": block["input"]}
                    })
                elif block["type"] == "tool_result":
                    ollama_messages.append({"role": "tool", "content": block["content"]})

            if text_parts or tool_calls:
                assistant_message = {"role": role, "content": "\n".join(text_parts)}
                if tool_calls:
                    assistant_message["tool_calls"] = tool_calls
                ollama_messages.append(assistant_message)

        return ollama_messages

    def send(self, messages: list, system_prompt: str, tools: list) -> dict:
        payload = {
            "model": self.model,
            "messages": self._to_ollama_messages(messages, system_prompt),
            "tools": self._to_ollama_tools(tools),
            "stream": False,
        }

        request = urllib.request.Request(
            f"{self.host}/api/chat",
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )

        try:
            with urllib.request.urlopen(request, timeout=300) as resp:
                body = json.loads(resp.read().decode("utf-8"))
        except urllib.error.URLError as e:
            print(
                f"Error: could not reach Ollama at {self.host} ({e}).\n"
                "Make sure the Ollama server is running (`ollama serve`) and the "
                f"model is pulled (`ollama pull {self.model}`)."
            )
            sys.exit(1)

        message = body.get("message", {})
        text = message.get("content", "")
        if text:
            print(text, end="", flush=True)

        content_blocks = []
        if text:
            content_blocks.append({"type": "text", "text": text})

        tool_calls = message.get("tool_calls") or []
        for call in tool_calls:
            function = call["function"]
            content_blocks.append({
                "type": "tool_use",
                "id": f"call_{uuid.uuid4().hex[:24]}",
                "name": function["name"],
                "input": function.get("arguments", {}),
            })

        stop_reason = "tool_use" if tool_calls else "end_turn"
        return {"content": content_blocks, "stop_reason": stop_reason}
