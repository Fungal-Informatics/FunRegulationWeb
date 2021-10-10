import { uuid } from "../helpers/utils";
import { User } from "../sdkgen/api-generated";
import { QUERY } from "../sql";

async function authUser(
	deviceId: string,
	email: string,
	password: string,
): Promise<void> {
	const user = await QUERY.users.getUserByEmailAndPassword(email, password);
	if (!user) throw Error("Email ou senha incorretos");
	await QUERY.users.setUserSessionKey(user.id, deviceId);
}

async function deauthUser(deviceId: string): Promise<void> {
	await QUERY.users.removeUserSessionKey(deviceId);
}

async function getUser(userId: string): Promise<User> {
	const user = await QUERY.users.getUser(userId);
	if (!user) throw new Error("Usuario nao existe");
	return user;
}

async function createUser(
	name: string,
	email: string,
	password: string,
): Promise<void> {
	const newId = uuid();

	// check if user exists
	const existingUser = await QUERY.users.getUserByEmail(email);
	if (existingUser) throw new Error("Já existe um usuario com esse email");

	// validate fields
	if (!email) throw Error("Email inválido");
	if (!password) throw Error("Senha inválida");

	// create user
	await QUERY.users.createUser({ id: newId, email, name, password });
}

export const users = {
	createUser,
	authUser,
	deauthUser,
	getUser,
};
